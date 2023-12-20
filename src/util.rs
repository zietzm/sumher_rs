use std::path::{Path, PathBuf};

use crate::io::tagging::TagInfo;

use anyhow::Result;
use csv::ReaderBuilder;
use polars::prelude::*;
use tokio::sync::Semaphore;

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
    let mut handles = Vec::new();

    for path in gwas_paths.iter().skip(1) {
        let sem = sem.clone();
        let first_series = first_series.clone();
        let path = path.clone();
        let handle = rt.spawn(async move {
            let permit = sem.acquire().await.unwrap();
            let series = read_predictors(path.as_path())?;
            drop(permit);
            Ok(series == *first_series)
        });
        handles.push(handle);
    }

    let mut aligned = true;
    for handle in handles {
        let result: Result<bool> = rt.block_on(handle)?;
        aligned = aligned && result?;
    }

    if aligned {
        let series = Series::new("Predictor", (*first_series).clone());
        Ok(Some(series.into_frame()))
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
