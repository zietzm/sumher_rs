use std::path::{Path, PathBuf};

use crate::io::tagging::TagInfo;

use anyhow::Result;
use polars::prelude::*;

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
fn read_predictors<P>(path: P) -> Result<Series>
where
    P: Into<PathBuf>,
{
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
pub fn check_predictors_aligned<P>(gwas_paths: &[P]) -> Result<Option<DataFrame>>
where
    P: Into<PathBuf> + Clone,
{
    let first_series = read_predictors(gwas_paths[0].clone())?;

    let aligned = gwas_paths
        .iter()
        .skip(1)
        .try_fold(true, |acc, x| -> Result<bool> {
            let series = read_predictors(x.clone())?;
            Ok(acc && (series == first_series))
        })?;

    if aligned {
        Ok(Some(first_series.into_frame()))
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
