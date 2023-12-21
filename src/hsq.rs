use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use crate::ffi::{solve_cors_wrapper, solve_sums_wrapper, SolveCorsResult, SolveSumsResult};
use crate::io::gwas::{sumstat_processor, sumstat_reader, AlignedGwasSumstats, RawGwasSumstats};
use crate::io::tagging::{read_tagfile, TagInfo};
use crate::util::{align_if_possible, check_predictors_aligned, RuntimeSetup};

use anyhow::Result;
use crossbeam_channel::Receiver;
use indicatif::ProgressBar;
use serde::{Serialize, Serializer};

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

pub fn compute_h2(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    output_root: &Path,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
    println!("Read tagfile");
    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;
    if !aligned {
        return Err(anyhow::anyhow!(
            "GWAS summary statistics files are not aligned!"
        ));
    }
    let tag_info = Arc::new(tag_info);
    let predictor_order = Arc::new(tag_info.predictor_order.clone());
    let gwas_paths = Arc::new(gwas_paths.to_vec());
    let output_root = Arc::new(output_root.to_path_buf());

    let pb = ProgressBar::new(gwas_paths.len() as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );
    let pb = Arc::new(Mutex::new(pb));

    let (raw_sender, raw_receiver) = crossbeam_channel::bounded::<RawGwasSumstats>(10);
    let (aligned_sender, aligned_receiver) = crossbeam_channel::bounded::<AlignedGwasSumstats>(10);

    let reader_process = std::thread::spawn(move || sumstat_reader(&gwas_paths, &raw_sender));
    println!("Started reader");

    let mut sumstat_workers = Vec::new();
    for _ in 0..runtime_setup.n_threads {
        let predictor_order = predictor_order.clone();
        let raw_receiver = raw_receiver.clone();
        let aligned_sender = aligned_sender.clone();
        sumstat_workers.push(std::thread::spawn(move || {
            sumstat_processor(&predictor_order, &raw_receiver, &aligned_sender)
        }));
    }
    println!("Started processors");

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
    println!("Started h2 workers");

    // Wait for workers to finish
    reader_process.join().unwrap()?;
    println!("Reader finished");

    for worker in sumstat_workers {
        worker.join().unwrap()?;
    }
    drop(aligned_sender);
    println!("Sumstat workers finished");

    for worker in ldak_workers {
        worker.join().unwrap()?;
    }
    println!("LDAK workers finished");

    Ok(())
}

fn h2_processor(
    sumstat_receiver: &Receiver<AlignedGwasSumstats>,
    tag_info: &TagInfo,
    output_root: &Path,
    progress: Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    for sumstats in sumstat_receiver {
        let output_stem = output_root
            .file_stem()
            .unwrap_or(OsStr::new("sumher_rs"))
            .to_str()
            .unwrap();

        let output_name = output_root
            .parent()
            .unwrap()
            .join(format!("{}.{}", output_stem, sumstats.phenotype))
            .to_str()
            .unwrap()
            .to_string();

        println!("Running h2 on {}", sumstats.phenotype);
        let progress_file = output_name.clone() + ".progress.txt";

        let result = solve_sums_wrapper(
            &tag_info.tag_vec,
            &sumstats.chisq,
            &sumstats.sample_sizes,
            &tag_info.cat_vec,
            &tag_info.category_info.ssums,
            &progress_file,
            None,
        )?;
        let partitions = format_heritability(&result, &tag_info.category_info.names);
        let output_path = output_name + ".hsq";
        write_results(&output_path, &partitions)?;
        println!("Wrote {}", sumstats.phenotype);
        progress.lock().unwrap().inc(1);
    }

    Ok(())
}

fn load_phenotypes_chunk(
    gwas_paths: &[PathBuf],
    tag_info: &Arc<TagInfo>,
    runtime_setup: &RuntimeSetup,
) -> Result<Vec<Arc<AlignedGwasSumstats>>> {
    let predictor_order = Arc::new(tag_info.predictor_order.clone());

    let mut sumstats = Vec::new();
    let mut sumstat_workers = Vec::new();

    let (raw_sender, raw_receiver) = crossbeam_channel::bounded::<RawGwasSumstats>(10);
    let (aligned_sender, aligned_receiver) = crossbeam_channel::bounded::<AlignedGwasSumstats>(10);

    let gwas_paths = Arc::new(gwas_paths.to_vec());
    let reader_process = std::thread::spawn(move || sumstat_reader(&gwas_paths, &raw_sender));

    for _ in 0..runtime_setup.n_threads {
        let predictor_order = predictor_order.clone();
        let raw_receiver = raw_receiver.clone();
        let aligned_sender = aligned_sender.clone();
        sumstat_workers.push(std::thread::spawn(move || {
            sumstat_processor(&predictor_order, &raw_receiver, &aligned_sender)
        }));
    }

    reader_process.join().unwrap()?;
    drop(aligned_sender);

    for worker in sumstat_workers {
        worker.join().unwrap()?;
    }

    while let Ok(sumstat) = aligned_receiver.recv() {
        sumstats.push(Arc::new(sumstat));
    }

    Ok(sumstats)
}

fn compute_rg_chunk(
    left_chunk: &[Arc<AlignedGwasSumstats>],
    right_chunk: &[Arc<AlignedGwasSumstats>],
    include_diag: bool,
    tag_info: &Arc<TagInfo>,
    output_root: &Path,
    progress: Arc<Mutex<ProgressBar>>,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let (sender, receiver) =
        crossbeam_channel::unbounded::<(Arc<AlignedGwasSumstats>, Arc<AlignedGwasSumstats>)>();

    // Put all combinations into the channel
    for (i, left) in left_chunk.iter().enumerate() {
        for (j, right) in right_chunk.iter().enumerate().skip(i) {
            if i == j && !include_diag {
                continue;
            }
            sender.send((left.clone(), right.clone()))?;
        }
    }
    drop(sender);

    let mut ldak_workers = Vec::new();
    for _ in 0..runtime_setup.n_threads {
        let tag_info = tag_info.clone();
        let receiver = receiver.clone();
        let progress = progress.clone();
        let output_root = output_root.to_path_buf();
        ldak_workers.push(std::thread::spawn(move || {
            rg_processor(&receiver, &tag_info, &output_root, &progress)
        }));
    }

    for worker in ldak_workers {
        worker.join().unwrap()?;
    }

    Ok(())
}

fn rg_processor(
    receiver: &Receiver<(Arc<AlignedGwasSumstats>, Arc<AlignedGwasSumstats>)>,
    tag_info: &Arc<TagInfo>,
    output_root: &Path,
    progress: &Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    for (left, right) in receiver {
        let output_stem = output_root
            .file_stem()
            .unwrap_or(OsStr::new("sumher_rs"))
            .to_str()
            .unwrap();

        let output_name = output_root
            .parent()
            .unwrap()
            .join(format!(
                "{}.{}.{}",
                output_stem, left.phenotype, right.phenotype
            ))
            .to_str()
            .unwrap()
            .to_string();

        println!("Running rg on {} - {}", left.phenotype, right.phenotype);
        let progress_file = output_name.clone() + ".progress.txt";

        let result = solve_cors_wrapper(
            &tag_info.tag_vec,
            &left,
            &right,
            &tag_info.cat_vec,
            &tag_info.category_info.ssums,
            &progress_file,
            None,
        )?;

        let partitions = format_genetic_correlation(&result);
        let output_path = output_name + ".rg";
        write_results(&output_path, &partitions)?;
        println!("Wrote {} - {}", left.phenotype, right.phenotype);
        progress.lock().unwrap().inc(1);
    }

    Ok(())
}

pub fn compute_rg(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    output_root: &Path,
    chunk_size: usize,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
    println!("Read tagfile");
    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;
    if !aligned {
        return Err(anyhow::anyhow!(
            "GWAS summary statistics files are not aligned!"
        ));
    }
    let tag_info = Arc::new(tag_info);

    let n_total = (gwas_paths.len() * (gwas_paths.len() - 1) / 2) as u64;
    let pb = ProgressBar::new(n_total);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );
    let pb = Arc::new(Mutex::new(pb));

    let chunks = gwas_paths.chunks(chunk_size / 2).collect::<Vec<_>>();

    for (i, left_chunk) in chunks.iter().enumerate() {
        let left_chunk = left_chunk.to_vec();
        println!("Loading chunk {:?}", left_chunk);
        let left_sumstats = load_phenotypes_chunk(&left_chunk, &tag_info, runtime_setup)?;

        compute_rg_chunk(
            &left_sumstats,
            &left_sumstats,
            false,
            &tag_info,
            output_root,
            pb.clone(),
            runtime_setup,
        )?;

        // do_rg(&left_sumstats, &left_sumstats, false)?;

        for (j, right_chunk) in chunks.iter().enumerate() {
            if i >= j {
                continue;
            }

            let right_chunk = right_chunk.to_vec();
            let right_sumstats = load_phenotypes_chunk(&right_chunk, &tag_info, runtime_setup)?;

            compute_rg_chunk(
                &left_sumstats,
                &right_sumstats,
                true,
                &tag_info,
                output_root,
                pb.clone(),
                runtime_setup,
            )?;
        }
    }

    Ok(())
}
