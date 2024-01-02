use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use crate::db::{DbConnection, Progress};
use crate::ffi::{solve_cors_wrapper, solve_sums_wrapper, SolveCorsResult, SolveSumsResult};
use crate::io::gwas::{sumstat_processor, sumstat_reader, AlignedGwasSumstats, RawGwasSumstats};
use crate::io::tagging::{read_tagfile, TagInfo};
use crate::util::{align_if_possible, check_predictors_aligned, make_progressbar, RuntimeSetup};

use anyhow::{anyhow, Result};
use crossbeam_channel::{Receiver, Sender};
use indicatif::ProgressBar;
use serde::{Serialize, Serializer};

fn round_serialize<S>(x: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_str(&format!("{:.6}", x))
}

#[derive(Debug, Serialize)]
pub struct HsqResult {
    pub phenotype: String,
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,
}

#[derive(Debug, Serialize)]
pub struct RgResult {
    pub phenotype1: String,
    pub phenotype2: String,
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,
}

fn format_heritability(
    phenotype: &str,
    result: &SolveSumsResult,
    category_names: &[String],
) -> Vec<HsqResult> {
    let n_categories = (result.stats.len() - 3) / 6;
    let mut partitions = Vec::new();

    category_names.iter().enumerate().for_each(|(i, name)| {
        let estimate = result.stats[i];
        let se = result.stats[2 * n_categories + 1 + i];

        partitions.push(HsqResult {
            phenotype: phenotype.to_string(),
            component: name.clone(),
            estimate,
            se,
        })
    });

    partitions.push(HsqResult {
        phenotype: phenotype.to_string(),
        component: "Her_All".to_string(),
        estimate: result.stats[n_categories],
        se: result.stats[3 * n_categories + 1],
    });

    partitions
}

fn format_genetic_correlation(
    phenotype1: String,
    phenotype2: String,
    result: &SolveCorsResult,
) -> Vec<RgResult> {
    let num_parts = result.num_parts as usize;
    let total = 2 * (result.gcon + result.cept) as usize + 3 * num_parts + 1;
    let total2 = num_parts + (result.gcon + result.cept) as usize;
    let se_idx = 2 * total + 4 + num_parts;

    let mut partitions = vec![
        RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Her1_All".to_string(),
            estimate: result.stats[total],
            se: result.stats[se_idx],
        },
        RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Her2_All".to_string(),
            estimate: result.stats[total + 1],
            se: result.stats[se_idx + 1],
        },
        RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Coher_All".to_string(),
            estimate: result.stats[total + 2],
            se: result.stats[se_idx + 2],
        },
        RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Cor_All".to_string(),
            estimate: result.stats[total + 3],
            se: result.stats[se_idx + 3],
        },
    ];

    if result.gcon == 1 {
        partitions.push(RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Scaling1".to_string(),
            estimate: result.stats[num_parts],
            se: result.stats[2 * num_parts + 4 + total],
        });
        partitions.push(RgResult {
            phenotype1: phenotype1.clone(),
            phenotype2: phenotype2.clone(),
            component: "Scaling2".to_string(),
            estimate: result.stats[total2 + num_parts],
            se: result.stats[2 * num_parts + 4 + total + total2],
        });
    }

    partitions.push(RgResult {
        phenotype1: phenotype1.clone(),
        phenotype2: phenotype2.clone(),
        component: "Overlap".to_string(),
        estimate: result.stats[2 * total2 + num_parts],
        se: result.stats[2 * num_parts + 4 + total + 2 * total2],
    });

    partitions
}

fn h2_result_writer(
    receiver: &Receiver<HsqResult>,
    progress: &Arc<Mutex<ProgressBar>>,
    conn: &Mutex<DbConnection>,
) -> Result<()> {
    for partition in receiver {
        conn.lock().unwrap().write_h2(&partition)?;
        progress.lock().unwrap().inc(1);
    }

    Ok(())
}

fn rg_result_writer(
    receiver: &Receiver<RgResult>,
    progress: &Arc<Mutex<ProgressBar>>,
    conn: &Mutex<DbConnection>,
) -> Result<()> {
    for partition in receiver {
        conn.lock().unwrap().write_rg(&partition)?;
        progress.lock().unwrap().inc(1);
    }

    Ok(())
}

fn remove_computed_h2(gwas_paths: &[PathBuf], computed: &Progress) -> Vec<PathBuf> {
    gwas_paths
        .iter()
        .filter(|x| !computed.finished_h2.contains(x.to_str().unwrap()))
        .cloned()
        .collect()
}

pub fn compute_h2(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let n_original = gwas_paths.len() as u64;
    let pb = make_progressbar(n_original);

    let gwas_paths = remove_computed_h2(gwas_paths, &runtime_setup.computed);
    pb.lock().unwrap().inc(n_original - gwas_paths.len() as u64);

    let alignment_info = check_predictors_aligned(&gwas_paths, runtime_setup.skip_alignment_check)?;
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;
    if !aligned {
        return Err(anyhow!("GWAS summary files are not aligned!"));
    }
    let tag_info = Arc::new(tag_info);
    let predictor_order = Arc::new(tag_info.predictor_order.clone());

    let n = runtime_setup.chunk_size / 2;
    let (raw_sender, raw_receiver) = crossbeam_channel::bounded::<RawGwasSumstats>(n);
    let (aligned_sender, aligned_receiver) = crossbeam_channel::bounded::<AlignedGwasSumstats>(n);
    let (h2_sender, h2_receiver) = crossbeam_channel::unbounded::<HsqResult>();

    let reader_process = std::thread::spawn(move || sumstat_reader(&gwas_paths, &raw_sender));
    let conn = runtime_setup.db_conn.clone();
    let writer_process = std::thread::spawn(move || h2_result_writer(&h2_receiver, &pb, &conn));

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
        let tag_info = tag_info.clone();
        let aligned_receiver = aligned_receiver.clone();
        let h2_sender = h2_sender.clone();
        ldak_workers.push(std::thread::spawn(move || {
            h2_processor(&aligned_receiver, &tag_info, &h2_sender)
        }));
    }

    // Wait for workers to finish
    reader_process.join().unwrap()?;

    for worker in sumstat_workers {
        worker.join().unwrap()?;
    }
    drop(aligned_sender);

    for worker in ldak_workers {
        worker.join().unwrap()?;
    }
    drop(h2_sender);

    writer_process.join().unwrap()?;

    Ok(())
}

fn h2_processor(
    sumstat_receiver: &Receiver<AlignedGwasSumstats>,
    tag_info: &TagInfo,
    result_sender: &Sender<HsqResult>,
) -> Result<()> {
    let names = &tag_info.category_info.names;

    for sumstats in sumstat_receiver {
        let result = solve_sums_wrapper(
            &tag_info.tag_vec,
            &sumstats.chisq,
            &sumstats.sample_sizes,
            &tag_info.cat_vec,
            &tag_info.category_info.ssums,
            None,
        )?;
        let partitions = format_heritability(&sumstats.phenotype, &result, names);
        for partition in partitions {
            result_sender.send(partition)?;
        }
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

    let (raw_sender, raw_receiver) = crossbeam_channel::unbounded::<RawGwasSumstats>();
    let (aligned_sender, aligned_receiver) = crossbeam_channel::unbounded::<AlignedGwasSumstats>();

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

    for worker in sumstat_workers {
        worker.join().unwrap()?;
    }
    drop(aligned_sender);

    while let Ok(sumstat) = aligned_receiver.recv() {
        sumstats.push(Arc::new(sumstat));
    }

    Ok(sumstats)
}

fn compute_rg_chunk(
    left_chunk: &[Arc<AlignedGwasSumstats>],
    right_chunk: &[Arc<AlignedGwasSumstats>],
    only_upper_tri: bool,
    result_sender: &Sender<RgResult>,
    tag_info: &Arc<TagInfo>,
    progress: &Arc<Mutex<ProgressBar>>,
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let (sender, receiver) =
        crossbeam_channel::unbounded::<(Arc<AlignedGwasSumstats>, Arc<AlignedGwasSumstats>)>();

    // Put all combinations into the channel
    for (i, left) in left_chunk.iter().enumerate() {
        for (j, right) in right_chunk.iter().enumerate() {
            if i >= j && only_upper_tri {
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
        let result_sender = result_sender.clone();
        let computed = runtime_setup.computed.clone();
        let pb = progress.clone();
        ldak_workers.push(std::thread::spawn(move || {
            rg_processor(&receiver, &tag_info, &result_sender, &computed, &pb)
        }));
    }

    for worker in ldak_workers {
        worker.join().unwrap()?;
    }

    Ok(())
}

fn rg_pre_computed(phenotype1: &str, phenotype2: &str, computed: &Progress) -> bool {
    computed
        .finished_rg
        .contains(&(phenotype1.to_string(), phenotype2.to_string()))
        || computed
            .finished_rg
            .contains(&(phenotype2.to_string(), phenotype1.to_string()))
}

fn rg_processor(
    receiver: &Receiver<(Arc<AlignedGwasSumstats>, Arc<AlignedGwasSumstats>)>,
    tag_info: &Arc<TagInfo>,
    result_sender: &Sender<RgResult>,
    precomputed: &Progress,
    progress: &Arc<Mutex<ProgressBar>>,
) -> Result<()> {
    for (left, right) in receiver {
        if rg_pre_computed(&left.phenotype, &right.phenotype, precomputed) {
            progress.lock().unwrap().inc(1);
            continue;
        }

        let result = solve_cors_wrapper(
            &tag_info.tag_vec,
            &left,
            &right,
            &tag_info.cat_vec,
            &tag_info.category_info.ssums,
            None,
        );
        match result {
            Ok(result) => {
                let partitions = format_genetic_correlation(
                    left.phenotype.clone(),
                    right.phenotype.clone(),
                    &result,
                );

                for partition in partitions {
                    result_sender.send(partition)?;
                }
            }
            Err(e) => {
                println!(
                    "Error computing rg for {} and {}: {}",
                    left.phenotype, right.phenotype, e
                );
            }
        }
    }

    Ok(())
}

pub fn compute_rg(
    tag_path: &Path,
    gwas_paths: &[PathBuf],
    runtime_setup: &RuntimeSetup,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;
    let alignment_info = check_predictors_aligned(gwas_paths, runtime_setup.skip_alignment_check)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;
    if !aligned {
        return Err(anyhow::anyhow!(
            "GWAS summary statistics files are not aligned!"
        ));
    }
    let tag_info = Arc::new(tag_info);

    let (result_sender, result_receiver) = crossbeam_channel::unbounded::<RgResult>();
    let n_paths = gwas_paths.len() as u64;
    let pb = make_progressbar(n_paths * (n_paths - 1) / 2);

    let writer_process = {
        let pb = pb.clone();
        let conn = runtime_setup.db_conn.clone();
        std::thread::spawn(move || rg_result_writer(&result_receiver, &pb, &conn))
    };

    let chunks = gwas_paths
        .chunks(runtime_setup.chunk_size / 2)
        .collect::<Vec<_>>();

    for (i, left_chunk) in chunks.iter().enumerate() {
        let left_chunk = left_chunk.to_vec();
        let left_sumstats = load_phenotypes_chunk(&left_chunk, &tag_info, runtime_setup)?;
        let pb = pb.clone();

        compute_rg_chunk(
            &left_sumstats,
            &left_sumstats,
            true,
            &result_sender,
            &tag_info,
            &pb,
            runtime_setup,
        )?;

        for (j, right_chunk) in chunks.iter().enumerate() {
            if i >= j {
                continue;
            }

            let right_chunk = right_chunk.to_vec();
            let right_sumstats = load_phenotypes_chunk(&right_chunk, &tag_info, runtime_setup)?;
            let pb = pb.clone();

            compute_rg_chunk(
                &left_sumstats,
                &right_sumstats,
                false,
                &result_sender,
                &tag_info,
                &pb,
                runtime_setup,
            )?;
        }
    }

    drop(result_sender);
    writer_process.join().unwrap()?;

    Ok(())
}
