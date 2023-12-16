use polars::prelude::*;

pub fn read_gwas_result(filename: &str) -> Result<DataFrame, Box<dyn std::error::Error>> {
    let df = CsvReader::from_path(filename)?
        .has_header(true)
        .with_separator(b'\t')
        .with_ignore_errors(true)
        .finish()?;

    Ok(df)
}
