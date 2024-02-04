use std::path::Path;
use std::path::PathBuf;

use anyhow::Result;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;

use crate::align::{align_gwas_taginfo, check_aligned};
use crate::db::Progress;
use crate::ffi::{solve_sums_wrapper, SolveSumsResult};
use crate::io::common::round_serialize;
use crate::io::gwas::GwasSumstats;
use crate::io::tagging::TagInfo;
use crate::util::{make_progressbar, RuntimeSetup};

pub fn compute_h2<P>(
    tag_path: P,
    gwas_paths: &[PathBuf],
    runtime_setup: &RuntimeSetup,
) -> Result<()>
where
    P: AsRef<Path>,
{
    let gwas_paths = remove_computed_h2(gwas_paths, &runtime_setup.computed);
    if gwas_paths.is_empty() {
        return Ok(());
    }
    let tag_info = TagInfo::from_file(tag_path)?;
    let pb = make_progressbar(gwas_paths.len() as u64);

    for chunk in gwas_paths
        .iter()
        .chunks(runtime_setup.chunk_size)
        .into_iter()
    {
        let results = chunk
            .collect::<Vec<_>>()
            .par_iter()
            .map(|x| {
                let mut gwas = GwasSumstats::from_file(x)?; // TODO: error handling after pb
                let res = wrap_h2(&mut tag_info.clone(), &mut gwas);
                pb.lock().unwrap().inc(1);
                res
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        runtime_setup.db_conn.lock().unwrap().write_h2(&results)?;
    }

    Ok(())
}

fn remove_computed_h2(gwas_paths: &[PathBuf], computed: &Progress) -> Vec<PathBuf> {
    gwas_paths
        .iter()
        .filter(|x| {
            !computed
                .finished_h2
                .contains(x.file_name().unwrap().to_str().unwrap())
        })
        .cloned()
        .collect()
}

fn wrap_h2(tag_info: &mut TagInfo, gwas_sumstats: &mut GwasSumstats) -> Result<Vec<HsqResult>> {
    if !check_aligned(gwas_sumstats, tag_info) {
        (*gwas_sumstats, *tag_info) = align_gwas_taginfo(gwas_sumstats, tag_info);
    }

    let result = solve_sums_wrapper(
        &tag_info.tagging,
        &gwas_sumstats.chisq,
        &gwas_sumstats.sample_size,
        &tag_info.annotations,
        &tag_info.ssums,
        None,
    )?;

    let h2_results = format_heritability(
        &gwas_sumstats.phenotype,
        &result,
        &tag_info.annotation_names,
    );

    Ok(h2_results)
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
