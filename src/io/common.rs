use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Result};
use serde::Serializer;

pub fn get_delimiter<P>(path: P) -> Result<u8>
where
    P: AsRef<Path>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let first_line = reader.lines().next().unwrap()?;

    if !first_line.starts_with("Predictor") {
        bail!("First line does not start with 'Predictor'");
    }

    let delimiter = first_line.as_bytes()[9];

    Ok(delimiter)
}

pub fn round_serialize<S>(x: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_str(&format!("{:.6}", x))
}
