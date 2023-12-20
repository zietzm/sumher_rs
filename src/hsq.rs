use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use crate::ffi::{solve_cors_wrapper, solve_sums_wrapper, SolveCorsResult, SolveSumsResult};
use crate::io::gwas::{
    read_gwas_aligned, read_gwas_result, sumstat_processor, sumstat_reader, AlignedGwasSumstats,
    GwasResultLine,
};
use crate::io::tagging::{read_tagfile, TagInfo};
use crate::util::{align_if_possible, check_predictors_aligned, RuntimeSetup};

use anyhow::Result;
use crossbeam_channel::Receiver;
use indicatif::ProgressBar;
use itertools::Itertools;
use polars::prelude::*;
use rayon::prelude::*;
use serde::{Serialize, Serializer};
use tokio::sync::Semaphore;
use tokio::task::spawn_blocking;

fn round_serialize<S>(x: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_str(&format!("{:.6}", x))
}

#[derive(Debug, Serialize)]
struct HeritabilityPartition {
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,

    #[serde(serialize_with = "round_serialize")]
    pub influence: f64,

    #[serde(serialize_with = "round_serialize")]
    pub influence_se: f64,
}

#[derive(Debug, Serialize)]
struct GeneticCorrelationPartition {
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,
}

fn format_heritability(
    result: &SolveSumsResult,
    category_names: &[String],
) -> Vec<HeritabilityPartition> {
    let n_categories = (result.stats.len() - 3) / 6;
    let mut partitions = Vec::new();

    category_names.iter().enumerate().for_each(|(i, name)| {
        let estimate = result.stats[i];
        let se = result.stats[2 * n_categories + 1 + i];

        partitions.push(HeritabilityPartition {
            component: name.clone(),
            estimate,
            se,
            influence: result.influs[i] * estimate,
            influence_se: result.influs[i].abs() * se,
        })
    });

    let mut total_influence = 0.0;
    let mut total_influence_se = 0.0;
    for i in 0..n_categories {
        total_influence += result.influs[i] * result.stats[i];
        for j in 0..n_categories {
            total_influence_se +=
                result.influs[i] * result.influs[j] * result.cohers[i + j * n_categories];
        }
    }

    partitions.push(HeritabilityPartition {
        component: "Her_All".to_string(),
        estimate: result.stats[n_categories],
        se: result.stats[3 * n_categories + 1],
        influence: total_influence,
        influence_se: total_influence_se.sqrt(),
    });

    partitions
}

fn format_genetic_correlation(result: &SolveCorsResult) -> Vec<GeneticCorrelationPartition> {
    let num_parts = result.num_parts as usize;
    let total = 2 * (result.gcon + result.cept) as usize + 3 * num_parts + 1;
    let total2 = num_parts + (result.gcon + result.cept) as usize;
    let se_idx = 2 * total + 4 + num_parts;

    let mut partitions = vec![
        GeneticCorrelationPartition {
            component: "Her1_All".to_string(),
            estimate: result.stats[total],
            se: result.stats[se_idx],
        },
        GeneticCorrelationPartition {
            component: "Her2_All".to_string(),
            estimate: result.stats[total + 1],
            se: result.stats[se_idx + 1],
        },
        GeneticCorrelationPartition {
            component: "Coher_All".to_string(),
            estimate: result.stats[total + 2],
            se: result.stats[se_idx + 2],
        },
        GeneticCorrelationPartition {
            component: "Cor_All".to_string(),
            estimate: result.stats[total + 3],
            se: result.stats[se_idx + 3],
        },
    ];

    if result.gcon == 1 {
        partitions.push(GeneticCorrelationPartition {
            component: "Scaling1".to_string(),
            estimate: result.stats[num_parts],
            se: result.stats[2 * num_parts + 4 + total],
        });
        partitions.push(GeneticCorrelationPartition {
            component: "Scaling2".to_string(),
            estimate: result.stats[total2 + num_parts],
            se: result.stats[2 * num_parts + 4 + total + total2],
        });
    }

    partitions.push(GeneticCorrelationPartition {
        component: "Overlap".to_string(),
        estimate: result.stats[2 * total2 + num_parts],
        se: result.stats[2 * num_parts + 4 + total + 2 * total2],
    });

    partitions
}

fn write_results<P, T: Serialize>(path: &P, results: &[T]) -> Result<()>
where
    P: AsRef<Path> + AsRef<std::ffi::OsStr> + ?Sized,
{
    // Write using csv and serde serialization
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;

    for result in results {
        writer.serialize(result)?;
    }

    writer.flush()?;
    Ok(())
}

// struct SumherInput {
//     tagging: Vec<f64>,
//     gwas_sumstats: AlignedGwasSumstats,
//     category_values: Vec<Vec<f64>>,
// }
//
// impl SumherInput {
//     fn new(
//         tagging: Vec<f64>,
//         gwas_sumstats: AlignedGwasSumstats,
//         category_values: Vec<Vec<f64>>,
//     ) -> Self {
//         Self {
//             tagging,
//             gwas_sumstats,
//             category_values,
//         }
//     }
//
//     fn from_aligned(tag_info: &TagInfo, gwas_df: &DataFrame) -> Result<Self> {
//         let tagging = tag_info.tag_vec.clone();
//         let category_values = tag_info.cat_vec.clone();
//         let gwas_sumstats = AlignedGwasSumstats::from_dataframe(gwas_df)?;
//         Ok(Self::new(tagging, gwas_sumstats, category_values))
//     }
//
//     fn from_misaligned(tag_info: &TagInfo, gwas_df: &DataFrame) -> Result<Self> {
//         let full_df = gwas_df
//             .join(
//                 &tag_info.df,
//                 ["Predictor"],
//                 ["Predictor"],
//                 JoinArgs::new(JoinType::Left),
//             )?
//             .lazy()
//             .with_column(col("Z").pow(2).alias("chisq"))
//             .collect()?;
//
//         let tagging = tag_info.tag_vec.clone();
//         let category_values = tag_info.cat_vec.clone();
//         let gwas_sumstats = AlignedGwasSumstats::from_dataframe(&full_df)?;
//         Ok(Self::new(tagging, gwas_sumstats, category_values))
//     }
//
//     fn from_gwas_tag_info(tag_info: &TagInfo, gwas_df: &DataFrame, aligned: bool) -> Result<Self> {
//         if aligned {
//             Self::from_aligned(tag_info, gwas_df)
//         } else {
//             Self::from_misaligned(tag_info, gwas_df)
//         }
//     }
// }

struct SumcorsInput {
    tagging: Vec<f64>,
    gwas_sumstats1: AlignedGwasSumstats,
    gwas_sumstats2: AlignedGwasSumstats,
    category_values: Vec<Vec<f64>>,
}

impl SumcorsInput {
    fn from_gwas_tag_info(
        tag_info: &TagInfo,
        gwas_df1: &DataFrame,
        gwas_df2: &DataFrame,
    ) -> Result<Self> {
        let full_df = tag_info
            .df
            .join(
                gwas_df1,
                ["Predictor"],
                ["Predictor"],
                JoinArgs::new(JoinType::Inner),
            )?
            .join(
                gwas_df2,
                ["Predictor"],
                ["Predictor"],
                JoinArgs::new(JoinType::Inner),
            )?
            .lazy()
            .with_columns([
                col("Z").pow(2).alias("chisq"),
                col("Z_right").pow(2).alias("chisq_right"),
            ])
            .collect()?;

        let gwas_df1_aligned = full_df.select(["Predictor", "n", "Z", "chisq"])?;
        let mut gwas_df2_aligned =
            full_df.select(["Predictor", "n_right", "Z_right", "chisq_right"])?;
        gwas_df2_aligned.set_column_names(&["Predictor", "n", "Z", "chisq"])?;

        let gwas_sumstats1 = AlignedGwasSumstats::from_dataframe(&gwas_df1_aligned)?;
        let gwas_sumstats2 = AlignedGwasSumstats::from_dataframe(&gwas_df2_aligned)?;

        let tag_info = TagInfo::from_dataframe(full_df, tag_info.category_info.clone())?;

        Ok(Self {
            tagging: tag_info.tag_vec,
            gwas_sumstats1,
            gwas_sumstats2,
            category_values: tag_info.cat_vec,
        })
    }
}

/// Compute heritability using LDAK across many files in parallel
// pub fn compute_hsq_parallel(
//     tag_path: &Path,
//     gwas_paths: &[PathBuf],
//     output_root: &Path,
//     runtime_setup: &RuntimeSetup,
// ) -> Result<()> {
//     let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
//
//     let pb = ProgressBar::new(gwas_paths.len() as u64);
//     pb.set_style(
//         indicatif::ProgressStyle::default_bar()
//             .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
//             .progress_chars("##-"),
//     );
//     let pb = Arc::new(Mutex::new(pb));
//     let output_root = Arc::new(output_root.to_path_buf());
//
//     let alignment_info = check_predictors_aligned(gwas_paths)?;
//     let aligned = align_if_possible(&mut tag_info, alignment_info)?;
//
//     let gwas_paths = gwas_paths
//         .iter()
//         .map(|x| Arc::new(x.clone()))
//         .collect::<Vec<Arc<PathBuf>>>();
//
//     if aligned {
//         compute_hsq_aligned(runtime_setup, &tag_info, &gwas_paths, output_root, pb)
//     } else {
//         compute_hsq_misaligned(runtime_setup, &tag_info, &gwas_paths, output_root, pb)
//     }
// }
//
// fn compute_hsq_aligned(
//     runtime_setup: &RuntimeSetup,
//     tag_info: &TagInfo,
//     gwas_paths: &[Arc<PathBuf>],
//     output_root: Arc<PathBuf>,
//     progress: Arc<Mutex<ProgressBar>>,
// ) -> Result<()> {
//     let tagging = Arc::new(tag_info.tag_vec.clone());
//     let category_values = Arc::new(tag_info.cat_vec.clone());
//     let category_contribs = Arc::new(tag_info.category_info.ssums.clone());
//     let category_names = Arc::new(tag_info.category_info.names.clone());
//
//     let tasks = gwas_paths
//         .par_iter()
//         .map(|x| {
//             let sem = runtime_setup.semaphore.clone();
//             let tag = tagging.clone();
//             let cat_val = category_values.clone();
//             let cat_con = category_contribs.clone();
//             let cat_names = category_names.clone();
//             let out = output_root.clone();
//             let x = x.clone();
//             let pb = progress.clone();
//             runtime_setup.runtime.spawn(async move {
//                 let result = h2_aligned(&tag, &cat_val, &cat_con, &cat_names, &x, &out, sem).await;
//                 pb.lock().unwrap().inc(1);
//                 result
//             })
//         })
//         .collect::<Vec<_>>();
//
//     for task in tasks {
//         let result = runtime_setup.runtime.block_on(task);
//         match result {
//             Ok(_) => {}
//             Err(e) => println!("Error: {}", e),
//         }
//     }
//
//     Ok(())
// }
//
// fn compute_hsq_misaligned(
//     runtime_setup: &RuntimeSetup,
//     tag_info: &TagInfo,
//     gwas_paths: &[Arc<PathBuf>],
//     output_root: Arc<PathBuf>,
//     progress: Arc<Mutex<ProgressBar>>,
// ) -> Result<()> {
//     let tag_info = Arc::new(tag_info.clone());
//
//     let tasks = gwas_paths
//         .par_iter()
//         .map(|x| {
//             let sem = runtime_setup.semaphore.clone();
//             let tag = tag_info.clone();
//             let out = output_root.clone();
//             let x = x.clone();
//             let pb = progress.clone();
//             runtime_setup.runtime.spawn(async move {
//                 let result = h2_misaligned(&tag, &x, &out, sem).await;
//                 pb.lock().unwrap().inc(1);
//                 result
//             })
//         })
//         .collect::<Vec<_>>();
//
//     for task in tasks {
//         let result = runtime_setup.runtime.block_on(task);
//         match result {
//             Ok(_) => {}
//             Err(e) => println!("Error: {}", e),
//         }
//     }
//
//     Ok(())
// }
//
/// Compute genetic correlations between all pairs of phenotypes using LDAK
pub fn compute_rg_parallel(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    output_root: &Path,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;

    let combinations = gwas_paths
        .iter()
        .map(|x| Arc::new(x.clone()))
        .tuple_combinations::<(_, _)>()
        .collect::<Vec<(Arc<PathBuf>, Arc<PathBuf>)>>();

    let pb = ProgressBar::new(combinations.len() as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );
    let pb = Arc::new(Mutex::new(pb));
    let output_root = Arc::new(output_root.to_path_buf());

    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;

    if aligned {
        compute_rg_aligned(runtime_setup, &tag_info, &combinations, output_root, pb)
    } else {
        compute_rg_misaligned(runtime_setup, &tag_info, &combinations, output_root, pb)
    }
}

fn compute_rg_misaligned(
    runtime_setup: &RuntimeSetup,
    tag_info: &TagInfo,
    combinations: &[(Arc<PathBuf>, Arc<PathBuf>)],
    output_root: Arc<PathBuf>,
    progress: Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    let tag_info = Arc::new(tag_info.clone());

    let tasks = combinations
        .par_iter()
        .map(|(x, y)| {
            let sem_clone = runtime_setup.semaphore.clone();
            let tag_info = tag_info.clone();
            let output_root = output_root.clone();
            let x = x.clone();
            let y = y.clone();
            let pb = progress.clone();
            runtime_setup.runtime.spawn(async move {
                let result = rg_misaligned(&tag_info, x, y, output_root.as_path(), sem_clone).await;
                pb.lock().unwrap().inc(1);
                result
            })
        })
        .collect::<Vec<_>>();

    for task in tasks {
        let result = runtime_setup.runtime.block_on(task);
        match result {
            Ok(_) => {}
            Err(e) => println!("Error: {}", e),
        }
    }

    Ok(())
}

fn compute_rg_aligned(
    runtime_setup: &RuntimeSetup,
    tag_info: &TagInfo,
    combinations: &[(Arc<PathBuf>, Arc<PathBuf>)],
    output_root: Arc<PathBuf>,
    progress: Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    let tagging = Arc::new(tag_info.tag_vec.clone());
    let category_values = Arc::new(tag_info.cat_vec.clone());
    let category_contribs = Arc::new(tag_info.category_info.ssums.clone());

    let tasks = combinations
        .par_iter()
        .map(|(x, y)| {
            let sem = runtime_setup.semaphore.clone();
            let tag = tagging.clone();
            let cat_val = category_values.clone();
            let cat_con = category_contribs.clone();
            let out = output_root.clone();
            let x = x.clone();
            let y = y.clone();
            let pb = progress.clone();
            runtime_setup.runtime.spawn(async move {
                let result = rg_aligned(&tag, &cat_val, &cat_con, x, y, out, sem).await;
                pb.lock().unwrap().inc(1);
                result
            })
        })
        .collect::<Vec<_>>();

    for task in tasks {
        let result = runtime_setup.runtime.block_on(task);
        match result {
            Ok(_) => {}
            Err(e) => println!("Error: {}", e),
        }
    }

    Ok(())
}

// async fn h2_misaligned(
//     tag_info: &TagInfo,
//     gwas_path: &Path,
//     output_path: &Path,
//     semaphore: Arc<Semaphore>,
// ) -> Result<()> {
//     let permit = semaphore.acquire().await?;
//     let gwas_df = read_gwas_result(gwas_path)?;
//     let input_data = SumherInput::from_gwas_tag_info(tag_info, &gwas_df, false)?;
//     let progress_path = output_path.with_extension("progress.txt");
//
//     let result = solve_sums_wrapper(
//         &input_data.tagging,
//         &input_data.gwas_sumstats.chisq,
//         &input_data.gwas_sumstats.sample_sizes,
//         &input_data.category_values,
//         &tag_info.category_info.ssums,
//         progress_path.to_str().unwrap(),
//         None,
//     );
//
//     let partitions = format_heritability(&result, &tag_info.category_info.names);
//     write_results(output_path, &partitions)?;
//     drop(permit);
//
//     Ok(())
// }
//
// async fn h2_aligned(
//     tagging: &[f64],
//     category_values: &[Vec<f64>],
//     category_contribs: &[Vec<f64>],
//     category_names: &[String],
//     gwas_path: &Path,
//     output_path: &Path,
//     semaphore: Arc<Semaphore>,
// ) -> Result<()> {
//     let permit = semaphore.acquire().await?;
//     let gwas_stats = read_gwas_aligned(gwas_path)?;
//     let progress_path = output_path.with_extension("progress.txt");
//
//     let result = solve_sums_wrapper(
//         tagging,
//         &gwas_stats.chisq,
//         &gwas_stats.sample_sizes,
//         category_values,
//         category_contribs,
//         progress_path.to_str().unwrap(),
//         None,
//     );
//
//     let partitions = format_heritability(&result, category_names);
//     write_results(output_path, &partitions)?;
//     drop(permit);
//
//     Ok(())
// }

fn format_rg_output_path(gwas_path_1: &Path, gwas_path_2: &Path, output_root: &Path) -> PathBuf {
    let f1 = gwas_path_1.file_stem().unwrap().to_str().unwrap();
    let f2 = gwas_path_2.file_stem().unwrap().to_str().unwrap();
    let combo_name = format!("{}.{}.rg", f1, f2);
    output_root.with_extension(combo_name)
}

async fn rg_misaligned(
    tag_info: &TagInfo,
    gwas_path_1: Arc<PathBuf>,
    gwas_path_2: Arc<PathBuf>,
    output_root: &Path,
    semaphore: Arc<Semaphore>,
) -> Result<()> {
    let permit = semaphore.acquire().await?;
    let path_1 = gwas_path_1.clone();
    let path_2 = gwas_path_2.clone();
    let gwas_df1 = tokio::task::spawn_blocking(move || read_gwas_result(path_1.as_path()));
    let gwas_df2 = tokio::task::spawn_blocking(move || read_gwas_result(path_2.as_path()));
    let gwas_df1 = gwas_df1.await??;
    let gwas_df2 = gwas_df2.await??;

    let input_data = SumcorsInput::from_gwas_tag_info(tag_info, &gwas_df1, &gwas_df2)?;

    let output_path = format_rg_output_path(&gwas_path_1, &gwas_path_2, output_root);
    let progress_path = output_path.with_extension("progress.txt");

    let result = solve_cors_wrapper(
        &input_data.tagging,
        &input_data.gwas_sumstats1,
        &input_data.gwas_sumstats2,
        &input_data.category_values,
        &tag_info.category_info.ssums,
        progress_path.to_str().unwrap(),
        None,
    );

    let partitions = format_genetic_correlation(&result);

    write_results(&output_path, &partitions)?;
    drop(permit);

    Ok(())
}

async fn rg_aligned(
    tagging: &[f64],
    category_values: &[Vec<f64>],
    category_contribs: &[Vec<f64>],
    gwas_path_1: Arc<PathBuf>,
    gwas_path_2: Arc<PathBuf>,
    output_root: Arc<PathBuf>,
    semaphore: Arc<Semaphore>,
) -> Result<()> {
    let permit = semaphore.acquire().await?;
    let path_1 = gwas_path_1.clone();
    let path_2 = gwas_path_2.clone();
    let gwas_1 = spawn_blocking(move || read_gwas_aligned(path_1.as_path()));
    let gwas_2 = spawn_blocking(move || read_gwas_aligned(path_2.as_path()));
    let gwas_1 = gwas_1.await??;
    let gwas_2 = gwas_2.await??;

    let output_path = format_rg_output_path(&gwas_path_1, &gwas_path_2, &output_root);
    let progress_path = output_path.with_extension("progress.txt");

    let result = solve_cors_wrapper(
        tagging,
        &gwas_1,
        &gwas_2,
        category_values,
        category_contribs,
        progress_path.to_str().unwrap(),
        None,
    );

    let partitions = format_genetic_correlation(&result);

    write_results(&output_path, &partitions)?;
    drop(permit);

    Ok(())
}

pub fn compute_h2(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    output_root: &Path,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;
    if !aligned {
        return Err(anyhow::anyhow!(
            "GWAS summary statistics files are not aligned!"
        ));
    }
    let tag_info = Arc::new(tag_info);
    let predictor_order = Arc::new(tag_info.predictor_order.clone());
    let gwas_paths = gwas_paths
        .iter()
        .map(|x| Arc::new(x.clone()))
        .collect::<Vec<Arc<PathBuf>>>();
    let output_root = Arc::new(output_root.to_path_buf());

    let pb = ProgressBar::new(gwas_paths.len() as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );
    let pb = Arc::new(Mutex::new(pb));

    let (raw_sender, raw_receiver) = crossbeam_channel::unbounded::<Vec<GwasResultLine>>();
    let (aligned_sender, aligned_receiver) = crossbeam_channel::unbounded::<AlignedGwasSumstats>();

    std::thread::spawn(move || sumstat_reader(&gwas_paths, &raw_sender));

    let mut sumstat_workers = Vec::new();
    for _ in 0..runtime_setup.n_threads {
        let predictor_order = predictor_order.clone();
        let raw_receiver = raw_receiver.clone();
        let aligned_sender = aligned_sender.clone();
        sumstat_workers.push(std::thread::spawn(move || {
            sumstat_processor(&predictor_order, &raw_receiver, &aligned_sender)
        }));
    }

    let mut ldak_workers = Vec::new();
    for _ in 0..runtime_setup.n_threads {
        let out = output_root.clone();
        let tag_info = tag_info.clone();
        let aligned_receiver = aligned_receiver.clone();
        let pb = pb.clone();
        ldak_workers.push(std::thread::spawn(move || {
            h2_processor(&aligned_receiver, &tag_info, &out, pb)
        }));
    }

    Ok(())
}

fn h2_processor(
    sumstat_receiver: &Receiver<AlignedGwasSumstats>,
    tag_info: &TagInfo,
    output_root: &Path,
    progress: Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    for sumstats in sumstat_receiver {
        let result = solve_sums_wrapper(
            &tag_info.tag_vec,
            &sumstats.chisq,
            &sumstats.sample_sizes,
            &tag_info.category_info.ssums,
            &tag_info.category_info.ssums,
            output_root.to_str().unwrap(),
            None,
        );
        let partitions = format_heritability(&result, &tag_info.category_info.names);
        write_results(output_root, &partitions)?;

        progress.lock().unwrap().inc(1);
    }

    Ok(())
}
