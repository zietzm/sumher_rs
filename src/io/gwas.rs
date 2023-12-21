use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use crossbeam_channel::{Receiver, Sender};
use csv::ReaderBuilder;
use itertools::izip;
use polars::prelude::*;

use serde::Deserialize;

use super::tagging::TagInfo;

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
// pub fn read_gwas_aligned<P>(filename: &P) -> Result<AlignedGwasSumstats>
// where
//     P: AsRef<Path> + AsRef<std::ffi::OsStr> + ?Sized,
// {
//     let mut n = Vec::new();
//     let mut z = Vec::new();
//
//     let mut csv_reader = ReaderBuilder::new()
//         .has_headers(true)
//         .delimiter(b'\t')
//         .from_path(filename)?;
//
//     let n_col = csv_reader.headers()?.iter().position(|x| x == "n").unwrap();
//     let z_col = csv_reader.headers()?.iter().position(|x| x == "Z").unwrap();
//
//     for result in csv_reader.records() {
//         let record = result?;
//         n.push(record.get(n_col).unwrap().parse::<f64>()?);
//         z.push(record.get(z_col).unwrap().parse::<f64>()?);
//     }
//
//     let chisq = z.iter().map(|x| x.powi(2)).collect::<Vec<f64>>();
//     let rho = izip!(chisq.iter(), n.iter(), z.iter())
//         .map(|(chi, ni, zi)| zi.signum() * (chi / (chi + ni)).sqrt())
//         .collect::<Vec<f64>>();
//
//     let phenotype: String = filename
//         .as_ref()
//         .file_stem()
//         .unwrap()
//         .to_str()
//         .unwrap()
//         .to_string();
//
//     Ok(AlignedGwasSumstats::new(phenotype, chisq, n, rho))
// }

pub struct AlignedGwasSumstats {
    pub phenotype: String,
    pub chisq: Vec<f64>,
    pub sample_sizes: Vec<f64>,
    pub rhos: Vec<f64>,
}

impl AlignedGwasSumstats {
    pub fn new(phenotype: String, chisq: Vec<f64>, sample_sizes: Vec<f64>, rhos: Vec<f64>) -> Self {
        Self {
            phenotype,
            chisq,
            sample_sizes,
            rhos,
        }
    }

    pub fn from_dataframe(df: &DataFrame, phenotype: String) -> Result<Self> {
        let z = df
            .column("Z")?
            .f64()?
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .context("Z column contains null values!")?;

        let chisq = z.iter().map(|x| x.powi(2)).collect::<Vec<_>>();

        let sample_sizes = df
            .column("n")?
            .cast(&DataType::Float64)?
            .f64()?
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .context("n column contains null values!")?;

        let rhos = izip!(z.iter(), chisq.iter(), sample_sizes.iter())
            .map(|(z, chisq, n)| z.signum() * (chisq / (chisq + n)).sqrt())
            .collect::<Vec<_>>();

        Ok(Self::new(phenotype, chisq, sample_sizes, rhos))
    }
}

#[derive(Clone, Debug, Deserialize)]
pub struct GwasResultLine {
    pub predictor: String,
    pub n: f64,
    pub z: f64,
}

pub struct RawGwasSumstats {
    pub phenotype: String,
    pub results: Vec<GwasResultLine>,
}

/// Read GWAS summary statistics from an LDAK-formatted file.
pub fn read_sumstats(path: &Path) -> Result<RawGwasSumstats> {
    let mut reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(path)?;

    let mut results = Vec::new();
    for result in reader.deserialize() {
        let record: GwasResultLine = result?;
        results.push(record);
    }

    let phenotype = path.file_stem().unwrap().to_str().unwrap().to_string();

    Ok(RawGwasSumstats { phenotype, results })
}

fn align_sumstats(
    predictor_order: &[String],
    raw_stats: &[GwasResultLine],
) -> Result<Vec<GwasResultLine>> {
    let aligned_stats = raw_stats
        .iter()
        .filter(|s| predictor_order.iter().any(|p| p == &s.predictor))
        .cloned()
        .collect::<Vec<GwasResultLine>>();

    Ok(aligned_stats)
}

pub fn process_sumstats(
    predictor_order: &[String],
    raw_stats: &RawGwasSumstats,
) -> Result<AlignedGwasSumstats> {
    // let is_aligned = predictor_order
    //     .iter()
    //     .zip(raw_stats.iter())
    //     .all(|(p, s)| p == &s.predictor);
    //
    // let aligned_stats = if is_aligned {
    //     raw_stats
    // } else {
    //     align_sumstats(predictor_order, raw_stats)?.as_slice()
    // };
    let aligned_stats = raw_stats; // For now let's just assume the stats are aligned

    let mut chisq = Vec::new();
    let mut sample_sizes = Vec::new();
    let mut rhos = Vec::new();

    for row in aligned_stats.results.iter() {
        let chisq_i = row.z.powi(2);
        chisq.push(chisq_i);
        sample_sizes.push(row.n);
        rhos.push(row.z.signum() * (chisq_i / (chisq_i + row.n)).sqrt());
    }

    println!("Processed {}", raw_stats.phenotype);

    Ok(AlignedGwasSumstats::new(
        raw_stats.phenotype.clone(),
        chisq,
        sample_sizes,
        rhos,
    ))
}

pub fn sumstat_reader(
    gwas_paths: &[Arc<PathBuf>],
    raw_channel: &Sender<RawGwasSumstats>,
) -> Result<()> {
    for path in gwas_paths {
        let results = read_sumstats(path)?;
        raw_channel.send(results)?;
    }
    Ok(())
}

pub fn sumstat_processor(
    predictor_order: &[String],
    raw_channel: &Receiver<RawGwasSumstats>,
    processed_channel: &Sender<AlignedGwasSumstats>,
) -> Result<()> {
    for raw_stats in raw_channel {
        let aligned_stats = process_sumstats(predictor_order, &raw_stats)?;
        processed_channel.send(aligned_stats)?;
    }
    Ok(())
}
