use std::path::Path;

use anyhow::Result;
use csv::ReaderBuilder;
use itertools::izip;
use polars::prelude::*;

use crate::hsq::AlignedGwasSumstats;

/// Read GWAS summary statistics from an LDAK-formatted file.
pub fn read_gwas_result<P>(filename: &P) -> Result<DataFrame>
where
    P: AsRef<Path> + AsRef<std::ffi::OsStr> + ?Sized,
{
    let df = CsvReader::from_path(filename)?
        .has_header(true)
        .with_separator(b'\t')
        .with_columns(Some(vec![
            "Predictor".to_string(),
            "n".to_string(),
            "Z".to_string(),
        ]))
        .with_ignore_errors(true)
        .finish()?;

    Ok(df)
}

/// Read GWAS summary statistics from an LDAK-formatted file when the predictors
/// are known to be aligned across files. This is faster than reading the file
/// into a DataFrame and joining the predictors.
pub fn read_gwas_aligned<P>(filename: &P) -> Result<AlignedGwasSumstats>
where
    P: AsRef<Path> + AsRef<std::ffi::OsStr> + ?Sized,
{
    let mut n = Vec::new();
    let mut z = Vec::new();

    let mut csv_reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(filename)?;

    let n_col = csv_reader.headers()?.iter().position(|x| x == "n").unwrap();
    let z_col = csv_reader.headers()?.iter().position(|x| x == "Z").unwrap();

    for result in csv_reader.records() {
        let record = result?;
        n.push(record.get(n_col).unwrap().parse::<f64>()?);
        z.push(record.get(z_col).unwrap().parse::<f64>()?);
    }

    let chisq = z.iter().map(|x| x.powi(2)).collect::<Vec<f64>>();
    let rho = izip!(chisq.iter(), n.iter(), z.iter())
        .map(|(chi, ni, zi)| zi.signum() * (chi / (chi + ni)).sqrt())
        .collect::<Vec<f64>>();

    Ok(AlignedGwasSumstats::new(chisq, n, rho))
}
