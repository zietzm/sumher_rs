use std::path::Path;

use anyhow::Result;
use polars::prelude::*;

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
