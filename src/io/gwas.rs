use std::fs::File;
use std::path::Path;

use anyhow::Result;
use csv::ReaderBuilder;
use itertools::izip;
use serde::Deserialize;

use crate::io::common::get_delimiter;

#[derive(Clone, Debug, Deserialize)]
pub struct GwasLine {
    #[serde(rename = "Predictor")]
    pub predictor: String,
    pub n: u32,
    #[serde(rename = "Z")]
    pub z: f64,
}

#[derive(Clone, Debug)]
pub struct GwasSumstats {
    pub phenotype: String,
    pub results: Vec<GwasLine>,

    pub predictors: Vec<String>,
    pub z_scores: Vec<f64>,
    pub chisq: Vec<f64>,
    pub sample_size: Vec<f64>,
    pub rhos: Vec<f64>,
}

impl GwasSumstats {
    pub fn new(phenotype: String, results: Vec<GwasLine>) -> Self {
        let predictors = results
            .iter()
            .map(|x| x.predictor.clone())
            .collect::<Vec<String>>();
        let z_scores = results.iter().map(|x| x.z).collect::<Vec<f64>>();
        let chisq = results.iter().map(|x| x.z * x.z).collect::<Vec<f64>>();
        let sample_size = results.iter().map(|x| x.n as f64).collect::<Vec<f64>>();

        let mut rhos = Vec::with_capacity(results.len());
        for (&z, &c, &n) in izip!(&z_scores, &chisq, &sample_size) {
            rhos.push(z.signum() * (c / (c + n)).sqrt());
        }

        GwasSumstats {
            phenotype,
            results,
            predictors,
            chisq,
            sample_size,
            z_scores,
            rhos,
        }
    }

    /// Read GWAS summary statistics from an LDAK-formatted file.
    pub fn from_file<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let delimiter = get_delimiter(&path)?;

        let file = File::open(&path)?;
        let mut reader = ReaderBuilder::new()
            .has_headers(true)
            .delimiter(delimiter)
            .from_reader(file);

        let mut results = reader
            .deserialize()
            .collect::<Result<Vec<GwasLine>, csv::Error>>()?;

        results.iter_mut().for_each(|x| {
            if x.z == 0.0 {
                x.z = 1e-6;
            }
        });

        let phenotype = path
            .as_ref()
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .to_string();

        Ok(GwasSumstats::new(phenotype, results))
    }
}
