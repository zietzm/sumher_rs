use anyhow::{anyhow, Result};
use csv::{ReaderBuilder, WriterBuilder};
use glob::glob;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

use crate::hsq::GeneticCorrelationPartition;
use crate::util::make_progressbar;

fn get_result_paths(root_glob: &Path) -> Result<Vec<PathBuf>> {
    let mut sumstats_paths = Vec::new();
    for path in glob(root_glob.to_str().unwrap())? {
        let path = path?;
        if path.is_dir() {
            for path in glob(&format!("{}/*.sumstats.gz", path.to_str().unwrap()))? {
                sumstats_paths.push(path?);
            }
        } else {
            sumstats_paths.push(path);
        }
    }
    Ok(sumstats_paths)
}

fn read_results(path: &PathBuf) -> Result<Vec<GeneticCorrelationPartition>> {
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let mut results = Vec::new();
    for result in reader.deserialize() {
        let result: GeneticCorrelationPartition = result?;
        results.push(result);
    }
    Ok(results)
}

fn phenotypes_from_path(path: &Path) -> Result<(String, String)> {
    let path = path.file_name().unwrap().to_str().unwrap();

    let re = Regex::new(r"(?:plink_white_british\.)(\w+)(?:\.glm\.linear)").unwrap();
    let mut phenotypes: Vec<String> = re.find_iter(path).map(|x| x.as_str().to_string()).collect();

    if phenotypes.len() != 2 {
        return Err(anyhow!("Could not find two phenotypes in path {}", path));
    }

    phenotypes.sort();

    let phenotype_1 = phenotypes[0].clone();
    let phenotype_2 = phenotypes[1].clone();
    Ok((phenotype_1, phenotype_2))
}

#[derive(Debug, Deserialize, Serialize)]
struct RgResultsRow {
    phenotype_1: String,
    phenotype_2: String,
    component: String,
    estimate: f64,
    se: f64,
}

pub fn collect_results(root: &Path, output_path: &Path) -> Result<()> {
    let mut results_paths = get_result_paths(root)?;
    println!("Found {} results files", results_paths.len());
    results_paths.sort();

    let pb = make_progressbar(results_paths.len() as u64);
    pb.lock().unwrap().set_message("Collecting results");

    let mut final_results = Vec::new();
    for path in results_paths {
        let (phenotype_1, phenotype_2) = phenotypes_from_path(&path)?;
        let results = read_results(&path)?;
        for result in results {
            final_results.push(RgResultsRow {
                phenotype_1: phenotype_1.clone(),
                phenotype_2: phenotype_2.clone(),
                component: result.component,
                estimate: result.estimate,
                se: result.se,
            });
        }
        pb.lock().unwrap().inc(1);
    }

    let mut writer = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(output_path)?;

    for result in final_results {
        writer.serialize(result)?;
    }

    Ok(())
}
