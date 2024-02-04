use std::path::{Path, PathBuf};

use anyhow::Result;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;

use crate::align::{align_gwas_gwas, check_aligned};
use crate::ffi::{solve_cors_wrapper, SolveCorsResult};
use crate::io::common::round_serialize;
use crate::io::gwas::GwasSumstats;
use crate::io::tagging::TagInfo;
use crate::util::{make_progressbar, RuntimeSetup};

pub fn compute_rg<P>(
    tag_path: P,
    gwas_paths: &[PathBuf],
    runtime_setup: &RuntimeSetup,
) -> Result<()>
where
    P: AsRef<Path>,
{
    let n_combos = gwas_paths.len() * (gwas_paths.len() - 1) / 2;
    let pb = make_progressbar(n_combos as u64);

    let tag_info = TagInfo::from_file(tag_path)?;

    let chunks: Vec<Vec<PathBuf>> = gwas_paths
        .iter()
        .chunks(runtime_setup.chunk_size / 2)
        .into_iter()
        .map(|x| x.cloned().collect::<Vec<_>>())
        .collect::<Vec<_>>();

    for (i, left_chunk) in chunks.iter().enumerate() {
        let left_stats = load_gwas_chunk(left_chunk)?;
        let left_results = left_stats
            .iter()
            .combinations(2)
            .collect::<Vec<_>>()
            .par_iter()
            .map(|x| wrap_rg(&mut tag_info.clone(), &mut x[0].clone(), &mut x[1].clone()))
            .collect::<Result<Vec<_>>>()?;

        runtime_setup
            .db_conn
            .lock()
            .unwrap()
            .write_rg(&left_results.into_iter().flatten().collect::<Vec<_>>())?;

        for right_chunk in chunks.iter().skip(i) {
            let right_stats = load_gwas_chunk(right_chunk)?;

            let combos = left_stats
                .iter()
                .cartesian_product(right_stats.iter())
                .filter(|(x, y)| {
                    !runtime_setup
                        .computed
                        .finished_rg
                        .contains(&(x.phenotype.clone(), y.phenotype.clone()))
                })
                .collect::<Vec<_>>();

            let results = combos
                .par_iter()
                .map(|(x, y)| {
                    let result =
                        wrap_rg(&mut tag_info.clone(), &mut (*x).clone(), &mut (*y).clone());
                    pb.lock().unwrap().inc(1);
                    result
                })
                .collect::<Result<Vec<_>>>()?;

            runtime_setup
                .db_conn
                .lock()
                .unwrap()
                .write_rg(&results.into_iter().flatten().collect::<Vec<_>>())?;
        }
    }

    Ok(())
}

fn load_gwas_chunk(files: &[PathBuf]) -> Result<Vec<GwasSumstats>> {
    files.iter().map(GwasSumstats::from_file).collect()
}

fn wrap_rg(
    tag_info: &mut TagInfo,
    gwas1: &mut GwasSumstats,
    gwas2: &mut GwasSumstats,
) -> Result<Vec<RgResult>> {
    if !check_aligned(gwas1, tag_info) || !check_aligned(gwas2, tag_info) {
        (*gwas1, *gwas2, *tag_info) = align_gwas_gwas(gwas1, gwas2, tag_info);
    }

    let result = solve_cors_wrapper(
        &tag_info.tagging,
        gwas1,
        gwas2,
        &tag_info.annotations,
        &tag_info.ssums,
        None,
    )?;

    Ok(format_genetic_correlation(
        gwas1.phenotype.clone(),
        gwas2.phenotype.clone(),
        &result,
    ))
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
