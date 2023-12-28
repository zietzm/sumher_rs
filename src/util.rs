use std::path::{Path, PathBuf};

use crate::io::tagging::TagInfo;

use anyhow::Result;
use csv::ReaderBuilder;
use indicatif::ProgressBar;
use polars::prelude::*;
use std::sync::{Arc, Mutex};

pub struct RuntimeSetup {
    pub n_threads: usize,
}

impl RuntimeSetup {
    pub fn new(n_threads: usize) -> Self {
        RuntimeSetup { n_threads }
    }
}

pub fn make_progressbar(n_total: u64) -> Arc<Mutex<ProgressBar>> {
    let pb = ProgressBar::new(n_total);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")
            .unwrap()
            .progress_chars("##-"),
    );
    Arc::new(Mutex::new(pb))
}

/// Reformat Plink summary statistics files for use with LDAK
pub fn format_plink_sumstats<P>(gwas_path: P, output_path: P) -> Result<()>
where
    P: Into<PathBuf> + AsRef<Path>,
{
    let mut df = CsvReader::from_path(gwas_path)?
        .has_header(true)
        .with_ignore_errors(true)
        .with_comment_char(Some(b'?'))
        .with_separator(b'\t')
        .finish()?;

    df = df
        .lazy()
        .select(&[
            col("ID").alias("Predictor"),
            col("A1"),
            col("OMITTED").alias("A2"),
            col("OBS_CT").alias("n"),
            col("T_STAT").alias("Z"),
        ])
        .collect()?;

    let mut file = std::fs::File::create(output_path)?;

    CsvWriter::new(&mut file)
        .include_header(true)
        .with_separator(b'\t')
        .finish(&mut df)
        .expect("Error writing output file");

    Ok(())
}

/// Read Predictors from LDAK-format GWAS summary statistics files
/// Predictors are the SNP IDs, either rsIDs or chr:pos
fn read_predictors(path: &Path) -> Result<Vec<String>> {
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let predictor_idx = reader
        .headers()?
        .iter()
        .position(|x| x == "Predictor")
        .ok_or_else(|| anyhow::anyhow!("Could not find Predictor column"))?;

    let mut predictors = Vec::new();
    for result in reader.records() {
        let record = result?;
        predictors.push(record.get(predictor_idx).unwrap().to_string());
    }
    Ok(predictors)
}

fn compare_predictors(path: &Path, comparison: &[String]) -> Result<bool> {
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let predictor_idx = reader
        .headers()?
        .iter()
        .position(|x| x == "Predictor")
        .ok_or_else(|| anyhow::anyhow!("Could not find Predictor column"))?;

    for (i, result) in reader.records().enumerate() {
        let record = result?;
        let predictor = record.get(predictor_idx).unwrap();
        if comparison[i] != predictor {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Check if predictors are aligned across GWAS summary statistics files
/// If so, return the shared predictors. Aligned predictors enable much
/// faster computation.
pub fn check_predictors_aligned(
    gwas_paths: &[PathBuf],
    skip_check: bool,
) -> Result<Option<DataFrame>> {
    let first_series = read_predictors(&gwas_paths[0])?;

    let pb = indicatif::ProgressBar::new((gwas_paths.len() - 1) as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );

    let mut aligned = true;
    if !skip_check {
        for path in gwas_paths.iter().skip(1) {
            let result = compare_predictors(path, &first_series)?;
            aligned = aligned && result;
            pb.inc(1);
        }
    }

    if aligned {
        let series = Series::new("Predictor", &first_series);
        Ok(Some(series.into_frame()))
    } else {
        Ok(None)
    }
}

/// Align predictors from the tag file with the predictors from the GWAS summary statistics files
pub fn align_if_possible(tag_info: &mut TagInfo, alignment: Option<DataFrame>) -> Result<bool> {
    match alignment {
        Some(shared_predictors) => {
            let new_df = tag_info
                .df
                .join(
                    &shared_predictors,
                    ["Predictor"],
                    ["Predictor"],
                    JoinArgs::new(JoinType::Inner),
                )
                .unwrap();
            tag_info.update_from_dataframe(new_df)?;
            Ok(true)
        }
        None => Ok(false),
    }
}
