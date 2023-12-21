use std::path::{Path, PathBuf};

use anyhow::Result;
use crossbeam_channel::{Receiver, Sender};
use csv::ReaderBuilder;

use serde::Deserialize;

#[derive(Clone, Debug, Deserialize)]
pub struct GwasResultLine {
    #[serde(rename = "Predictor")]
    pub predictor: String,
    pub n: f64,
    #[serde(rename = "Z")]
    pub z: f64,
}

pub struct RawGwasSumstats {
    pub phenotype: String,
    pub results: Vec<GwasResultLine>,
}

pub struct AlignedGwasSumstats {
    pub phenotype: String,
    pub chisq: Vec<f64>,
    pub sample_sizes: Vec<f64>,
    pub rhos: Vec<f64>,
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
    let is_aligned = predictor_order
        .iter()
        .zip(raw_stats.results.iter())
        .all(|(p, s)| p == &s.predictor);

    let aligned_vec: Vec<GwasResultLine>;
    let aligned_stats: &[GwasResultLine];

    if is_aligned {
        aligned_stats = &raw_stats.results;
    } else {
        println!("WARNING: Aligning {}", raw_stats.phenotype);
        aligned_vec = align_sumstats(predictor_order, &raw_stats.results)?;
        aligned_stats = &aligned_vec;
    }

    let mut chisq = Vec::new();
    let mut sample_sizes = Vec::new();
    let mut rhos = Vec::new();

    for row in aligned_stats.iter() {
        let chisq_i = row.z.powi(2);
        chisq.push(chisq_i);
        sample_sizes.push(row.n);
        rhos.push(row.z.signum() * (chisq_i / (chisq_i + row.n)).sqrt());
    }

    println!("Processed {}", raw_stats.phenotype);

    Ok(AlignedGwasSumstats {
        phenotype: raw_stats.phenotype.clone(),
        chisq,
        sample_sizes,
        rhos,
    })
}

pub fn sumstat_reader(gwas_paths: &[PathBuf], raw_channel: &Sender<RawGwasSumstats>) -> Result<()> {
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
