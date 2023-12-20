use std::path::{Path, PathBuf};

use crate::io::tagging::TagInfo;

use anyhow::{Context, Result};
use polars::prelude::*;
use tokio::sync::Semaphore;
use tokio::task::JoinSet;

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
fn read_predictors(path: &Path) -> Result<Series> {
    let df = CsvReader::from_path(path)?
        .has_header(true)
        .with_separator(b'\t')
        .with_columns(Some(vec!["Predictor".to_string()]))
        .finish()?;

    Ok(df.column("Predictor").unwrap().clone())
}

/// Check if predictors are aligned across GWAS summary statistics files
/// If so, return the shared predictors. Aligned predictors enable much
/// faster computation.
pub fn check_predictors_aligned(gwas_paths: &[PathBuf]) -> Result<Option<DataFrame>> {
    let first_series = read_predictors(&gwas_paths[0])?;
    let first_series = Arc::new(first_series);

    let gwas_paths = gwas_paths
        .iter()
        .map(|x| Arc::new(x.clone()))
        .collect::<Vec<Arc<PathBuf>>>();

    let sem = Arc::new(Semaphore::new(num_cpus::get()));
    let rt = tokio::runtime::Runtime::new()?;

    let mut set = JoinSet::new();
    for path in gwas_paths.iter().skip(1) {
        let sem = sem.clone();
        let first_series = first_series.clone();
        let path = path.clone();
        set.spawn(async move {
            let permit = sem.acquire().await.unwrap();
            let series = read_predictors(path.as_path())?;
            drop(permit);
            Ok(series == *first_series)
        });
    }

    let mut aligned = true;
    while let Some(result) = rt.block_on(set.join_next()) {
        let result: Result<bool> = result.context("Error reading predictors")?;
        aligned = aligned && result?;
    }

    if aligned {
        Ok(Some((*first_series).clone().into_frame()))
    } else {
        Ok(None)
    }
}

/// Align predictors from the tag file with the predictors from the GWAS summary statistics files
pub fn align_if_possible(tag_info: &mut TagInfo, alignment: Option<DataFrame>) -> Result<bool> {
    match alignment {
        Some(shared_predictors) => {
            println!("Predictors are aligned.");
            tag_info.df = tag_info
                .df
                .join(
                    &shared_predictors,
                    ["Predictor"],
                    ["Predictor"],
                    JoinArgs::new(JoinType::Inner),
                )
                .unwrap();
            Ok(true)
        }
        None => Ok(false),
    }
}
