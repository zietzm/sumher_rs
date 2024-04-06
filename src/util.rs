use indicatif::ProgressBar;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use anyhow::Result;
use csv::{ReaderBuilder, WriterBuilder};
use serde::{Deserialize, Serialize};

use crate::db::{DbConnection, Progress};

pub struct RuntimeSetup {
    pub n_threads: usize,
    pub chunk_size: usize,
    pub db_conn: Arc<Mutex<DbConnection>>,
    pub skip_alignment_check: bool,
    pub computed: Arc<Progress>,
}

impl RuntimeSetup {
    pub fn new(
        n_threads: usize,
        chunk_size: usize,
        db_conn: DbConnection,
        skip_alignment_check: bool,
        computed: Progress,
    ) -> Self {
        RuntimeSetup {
            n_threads,
            chunk_size,
            db_conn: Arc::new(Mutex::new(db_conn)),
            skip_alignment_check,
            computed: Arc::new(computed),
        }
    }
}

pub fn make_progressbar(n_total: u64) -> Arc<Mutex<ProgressBar>> {
    let pb = ProgressBar::new(n_total);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40} {pos:>7}/{len:7} ({eta}) {msg}")
            .unwrap()
            .progress_chars("##-"),
    );
    Arc::new(Mutex::new(pb))
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PlinkSumstats {
    #[serde(rename(serialize = "Predictor", deserialize = "ID"))]
    predictor: String,
    #[serde(rename = "A1")]
    a1: String,
    #[serde(rename(serialize = "A2", deserialize = "OMITTED"))]
    a2: String,
    #[serde(rename(serialize = "n", deserialize = "OBS_CT"))]
    n: i64,
    #[serde(rename(serialize = "Z", deserialize = "T_STAT"))]
    z: f64,
}

/// Reformat Plink summary statistics files for use with LDAK
pub fn format_plink_sumstats<P>(gwas_path: P, output_path: P) -> Result<()>
where
    P: Into<PathBuf> + AsRef<Path>,
{
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .buffer_capacity(8 * (1 << 20))
        .from_path(&gwas_path)?;

    let mut writer = WriterBuilder::new()
        .delimiter(b' ')
        .buffer_capacity(8 * (1 << 20))
        .from_path(&output_path)?;

    for record in reader.deserialize::<PlinkSumstats>() {
        writer.serialize(record?)?;
    }

    Ok(())
}
