use anyhow::Result;
use itertools::Itertools;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use tokio::runtime::Runtime;
use tokio::sync::Semaphore;

use crate::hsq::estimate_genetic_correlation_async;
use crate::hsq::estimate_heritability;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle};
use rayon::prelude::*;

use crate::io::tagging::read_tagfile;
use crate::io::tagging::TagInfo;

use polars::prelude::*;

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

fn align_if_possible(tag_info: &mut TagInfo, alignment: Option<DataFrame>) -> Result<bool> {
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

pub fn compute_hsq_parallel(
    tag_path: &Path,
    gwas_paths: &[&PathBuf],
    output_root: &Path,
    chunk_size: usize,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;

    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;

    gwas_paths
        .chunks(chunk_size)
        .progress_count((gwas_paths.len() as f32 / chunk_size as f32).ceil() as u64)
        .for_each(|batch| {
            batch.par_iter().for_each(|path| {
                let output_path = output_root
                    .clone()
                    .join(path.file_stem().unwrap())
                    .with_extension("hers");

                let result = estimate_heritability(&tag_info, path, &output_path, aligned);
                match result {
                    Ok(_) => {}
                    Err(e) => println!("Error: {}", e),
                }
            })
        });

    Ok(())
}

pub fn compute_rg_parallel(
    tag_path: &Path,
    gwas_paths: &[&PathBuf],
    output_root: &Path,
    n_permits: usize,
) -> Result<()> {
    let mut tag_info = read_tagfile(tag_path.to_str().unwrap())?;

    let alignment_info = check_predictors_aligned(gwas_paths)?;
    let aligned = align_if_possible(&mut tag_info, alignment_info)?;

    let combinations = gwas_paths
        .iter()
        .tuple_combinations::<(_, _)>()
        .map(|(x, y)| (*x, *y))
        .map(|(x, y)| (x.to_path_buf(), y.to_path_buf()))
        .collect::<Vec<_>>();

    let rt = Runtime::new()?;
    let sem = Arc::new(Semaphore::new(n_permits));
    let tag_info = Arc::new(tag_info.clone());
    let output_root = output_root.to_path_buf();

    let pb = ProgressBar::new(combinations.len() as u64);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} ({eta}) {msg}")?
            .progress_chars("##-"),
    );

    let pb = Arc::new(Mutex::new(pb));

    let tasks = combinations
        .par_iter()
        .map(|(x, y)| {
            let sem_clone = sem.clone();
            let tag_info = tag_info.clone();
            let output_root = output_root.clone();
            let x = x.clone();
            let y = y.clone();
            let pb = pb.clone();
            rt.spawn(async move {
                let result = estimate_genetic_correlation_async(
                    tag_info,
                    x,
                    y,
                    output_root,
                    aligned,
                    sem_clone,
                )
                .await;
                pb.lock().unwrap().inc(1);
                result
            })
        })
        .collect::<Vec<_>>();

    for task in tasks {
        let result = rt.block_on(task);
        match result {
            Ok(_) => {}
            Err(e) => println!("Error: {}", e),
        }
    }

    Ok(())
}
