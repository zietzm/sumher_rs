use std::path::{Path, PathBuf};

use crate::ffi::solve_cors_wrapper;
use crate::ffi::solve_sums_wrapper;
use crate::ffi::{SolveCorsResult, SolveSumsResult};
use crate::io::gwas::{read_gwas_aligned, read_gwas_result};
use crate::io::tagging::TagInfo;

use anyhow::Context;
use anyhow::Result;
use itertools::izip;
use ndarray::Axis;
use polars::prelude::*;
use serde::Serialize;
use serde::Serializer;
use tokio::sync::Semaphore;
use tokio::task::spawn_blocking;

fn round_serialize<S>(x: &f64, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_str(&format!("{:.6}", x))
}

#[derive(Debug, Serialize)]
pub struct HeritabilityPartition {
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,

    #[serde(serialize_with = "round_serialize")]
    pub influence: f64,

    #[serde(serialize_with = "round_serialize")]
    pub influence_se: f64,
}

#[derive(Debug, Serialize)]
pub struct GeneticCorrelationPartition {
    pub component: String,

    #[serde(serialize_with = "round_serialize")]
    pub estimate: f64,

    #[serde(serialize_with = "round_serialize")]
    pub se: f64,
}

fn format_heritability(
    result: &SolveSumsResult,
    category_names: &[String],
) -> Vec<HeritabilityPartition> {
    let n_categories = (result.stats.len() - 3) / 6;
    let mut partitions = Vec::new();

    category_names.iter().enumerate().for_each(|(i, name)| {
        let estimate = result.stats[i];
        let se = result.stats[2 * n_categories + 1 + i];

        partitions.push(HeritabilityPartition {
            component: name.clone(),
            estimate,
            se,
            influence: result.influs[i] * estimate,
            influence_se: result.influs[i].abs() * se,
        })
    });

    let mut total_influence = 0.0;
    let mut total_influence_se = 0.0;
    for i in 0..n_categories {
        total_influence += result.influs[i] * result.stats[i];
        for j in 0..n_categories {
            total_influence_se +=
                result.influs[i] * result.influs[j] * result.cohers[i + j * n_categories];
        }
    }

    partitions.push(HeritabilityPartition {
        component: "Her_All".to_string(),
        estimate: result.stats[n_categories],
        se: result.stats[3 * n_categories + 1],
        influence: total_influence,
        influence_se: total_influence_se.sqrt(),
    });

    partitions
}

fn format_genetic_correlation(result: &SolveCorsResult) -> Vec<GeneticCorrelationPartition> {
    let num_parts = result.num_parts as usize;
    let total = 2 * (result.gcon + result.cept) as usize + 3 * num_parts + 1;
    let total2 = num_parts + (result.gcon + result.cept) as usize;
    let se_idx = 2 * total + 4 + num_parts;

    let mut partitions = vec![
        GeneticCorrelationPartition {
            component: "Her1_All".to_string(),
            estimate: result.stats[total],
            se: result.stats[se_idx],
        },
        GeneticCorrelationPartition {
            component: "Her2_All".to_string(),
            estimate: result.stats[total + 1],
            se: result.stats[se_idx + 1],
        },
        GeneticCorrelationPartition {
            component: "Coher_All".to_string(),
            estimate: result.stats[total + 2],
            se: result.stats[se_idx + 2],
        },
        GeneticCorrelationPartition {
            component: "Cor_All".to_string(),
            estimate: result.stats[total + 3],
            se: result.stats[se_idx + 3],
        },
    ];

    if result.gcon == 1 {
        partitions.push(GeneticCorrelationPartition {
            component: "Scaling1".to_string(),
            estimate: result.stats[num_parts],
            se: result.stats[2 * num_parts + 4 + total],
        });
        partitions.push(GeneticCorrelationPartition {
            component: "Scaling2".to_string(),
            estimate: result.stats[total2 + num_parts],
            se: result.stats[2 * num_parts + 4 + total + total2],
        });
    }

    partitions.push(GeneticCorrelationPartition {
        component: "Overlap".to_string(),
        estimate: result.stats[2 * total2 + num_parts],
        se: result.stats[2 * num_parts + 4 + total + 2 * total2],
    });

    partitions
}

fn write_results<P, T: Serialize>(path: &P, results: &[T]) -> Result<()>
where
    P: AsRef<Path> + AsRef<std::ffi::OsStr> + ?Sized,
{
    // Write using csv and serde serialization
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').from_path(path)?;

    for result in results {
        writer.serialize(result)?;
    }

    writer.flush()?;
    Ok(())
}

pub struct AlignedGwasSumstats {
    pub chisq: Vec<f64>,
    pub sample_sizes: Vec<f64>,
    pub rhos: Vec<f64>,
}

impl AlignedGwasSumstats {
    pub fn new(chisq: Vec<f64>, sample_sizes: Vec<f64>, rhos: Vec<f64>) -> Self {
        Self {
            chisq,
            sample_sizes,
            rhos,
        }
    }

    pub fn from_dataframe(df: &DataFrame) -> Result<Self> {
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

        Ok(Self::new(chisq, sample_sizes, rhos))
    }
}

pub fn get_tagging_vec(df: &DataFrame) -> Result<Vec<f64>> {
    df.column("Tagging")?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .context("Tagging column contains null values!")
}

pub fn get_category_values_vec(df: &DataFrame, category_names: &[String]) -> Result<Vec<Vec<f64>>> {
    Ok(df
        .select(category_names.iter())?
        .to_ndarray::<Float64Type>(IndexOrder::Fortran)?
        .axis_iter(Axis(1))
        .map(|x| x.to_vec())
        .collect::<Vec<Vec<f64>>>())
}

struct SumherInput {
    tagging: Vec<f64>,
    gwas_sumstats: AlignedGwasSumstats,
    category_values: Vec<Vec<f64>>,
}

impl SumherInput {
    fn new(
        tagging: Vec<f64>,
        gwas_sumstats: AlignedGwasSumstats,
        category_values: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            tagging,
            gwas_sumstats,
            category_values,
        }
    }

    fn from_aligned(tag_info: &TagInfo, gwas_df: &DataFrame) -> Result<Self> {
        let tagging = get_tagging_vec(&tag_info.df)?;
        let category_values = get_category_values_vec(&tag_info.df, &tag_info.category_info.names)?;
        let gwas_sumstats = AlignedGwasSumstats::from_dataframe(gwas_df)?;
        Ok(Self::new(tagging, gwas_sumstats, category_values))
    }

    fn from_misaligned(tag_info: &TagInfo, gwas_df: &DataFrame) -> Result<Self> {
        let full_df = gwas_df
            .join(
                &tag_info.df,
                ["Predictor"],
                ["Predictor"],
                JoinArgs::new(JoinType::Left),
            )?
            .lazy()
            .with_column(col("Z").pow(2).alias("chisq"))
            .collect()?;

        let tagging = get_tagging_vec(&full_df)?;
        let category_values = get_category_values_vec(&full_df, &tag_info.category_info.names)?;
        let gwas_sumstats = AlignedGwasSumstats::from_dataframe(&full_df)?;
        Ok(Self::new(tagging, gwas_sumstats, category_values))
    }

    fn from_gwas_tag_info(tag_info: &TagInfo, gwas_df: &DataFrame, aligned: bool) -> Result<Self> {
        if aligned {
            Self::from_aligned(tag_info, gwas_df)
        } else {
            Self::from_misaligned(tag_info, gwas_df)
        }
    }
}

pub struct SumcorsInput {
    tagging: Vec<f64>,
    gwas_sumstats1: AlignedGwasSumstats,
    gwas_sumstats2: AlignedGwasSumstats,
    category_values: Vec<Vec<f64>>,
}

impl SumcorsInput {
    fn new(
        tagging: Vec<f64>,
        gwas_sumstats1: AlignedGwasSumstats,
        gwas_sumstats2: AlignedGwasSumstats,
        category_values: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            tagging,
            gwas_sumstats1,
            gwas_sumstats2,
            category_values,
        }
    }

    // TODO: This function can probably be removed?
    fn from_aligned(
        tag_info: &TagInfo,
        gwas_df1: &DataFrame,
        gwas_df2: &DataFrame,
    ) -> Result<Self> {
        let tagging = get_tagging_vec(&tag_info.df)?;
        let category_values = get_category_values_vec(&tag_info.df, &tag_info.category_info.names)?;

        let gwas_sumstats1 = AlignedGwasSumstats::from_dataframe(gwas_df1)?;
        let gwas_sumstats2 = AlignedGwasSumstats::from_dataframe(gwas_df2)?;

        Ok(Self::new(
            tagging,
            gwas_sumstats1,
            gwas_sumstats2,
            category_values,
        ))
    }

    fn from_misaligned(
        tag_info: &TagInfo,
        gwas_df1: &DataFrame,
        gwas_df2: &DataFrame,
    ) -> Result<Self> {
        let full_df = tag_info
            .df
            .join(
                gwas_df1,
                ["Predictor"],
                ["Predictor"],
                JoinArgs::new(JoinType::Inner),
            )?
            .join(
                gwas_df2,
                ["Predictor"],
                ["Predictor"],
                JoinArgs::new(JoinType::Inner),
            )?
            .lazy()
            .with_columns([
                col("Z").pow(2).alias("chisq"),
                col("Z_right").pow(2).alias("chisq_right"),
            ])
            .collect()?;

        let gwas_df1_aligned = full_df.select(["Predictor", "n", "Z", "chisq"])?;
        let mut gwas_df2_aligned =
            full_df.select(["Predictor", "n_right", "Z_right", "chisq_right"])?;
        gwas_df2_aligned.set_column_names(&["Predictor", "n", "Z", "chisq"])?;

        let gwas_sumstats1 = AlignedGwasSumstats::from_dataframe(&gwas_df1_aligned)?;
        let gwas_sumstats2 = AlignedGwasSumstats::from_dataframe(&gwas_df2_aligned)?;

        let tagging = get_tagging_vec(&full_df)?;
        let category_values = get_category_values_vec(&full_df, &tag_info.category_info.names)?;

        Ok(Self::new(
            tagging,
            gwas_sumstats1,
            gwas_sumstats2,
            category_values,
        ))
    }

    fn from_gwas_tag_info(
        tag_info: &TagInfo,
        gwas_df1: &DataFrame,
        gwas_df2: &DataFrame,
        aligned: bool,
    ) -> Result<Self> {
        if aligned {
            Self::from_aligned(tag_info, gwas_df1, gwas_df2)
        } else {
            Self::from_misaligned(tag_info, gwas_df1, gwas_df2)
        }
    }
}

pub fn estimate_heritability(
    tag_info: &TagInfo,
    gwas_path: &Path,
    output_path: &Path,
    aligned: bool,
) -> Result<()> {
    let gwas_df = read_gwas_result(gwas_path)?;

    let input_data = SumherInput::from_gwas_tag_info(tag_info, &gwas_df, aligned)?;
    let progress_path = output_path.with_extension("progress.txt");

    let result = solve_sums_wrapper(
        &input_data.tagging,
        &input_data.gwas_sumstats.chisq,
        &input_data.gwas_sumstats.sample_sizes,
        &input_data.category_values,
        &tag_info.category_info.ssums,
        progress_path.to_str().unwrap(),
        None,
    );

    let partitions = format_heritability(&result, &tag_info.category_info.names);
    write_results(output_path, &partitions)?;

    Ok(())
}

fn format_rg_output_path(gwas_path_1: &Path, gwas_path_2: &Path, output_root: &Path) -> PathBuf {
    let f1 = gwas_path_1.file_stem().unwrap().to_str().unwrap();
    let f2 = gwas_path_2.file_stem().unwrap().to_str().unwrap();
    let combo_name = format!("{}.{}.rg", f1, f2);
    output_root.with_extension(combo_name)
}

pub async fn rg_misaligned(
    tag_info: Arc<TagInfo>,
    gwas_path_1: Arc<PathBuf>,
    gwas_path_2: Arc<PathBuf>,
    output_root: Arc<PathBuf>,
    semaphore: Arc<Semaphore>,
) -> Result<()> {
    let permit = semaphore.acquire().await?;
    let path_1 = gwas_path_1.clone();
    let path_2 = gwas_path_2.clone();
    let gwas_df1 = tokio::task::spawn_blocking(move || read_gwas_result(path_1.as_path()));
    let gwas_df2 = tokio::task::spawn_blocking(move || read_gwas_result(path_2.as_path()));
    let gwas_df1 = gwas_df1.await??;
    let gwas_df2 = gwas_df2.await??;

    let input_data = SumcorsInput::from_gwas_tag_info(&tag_info, &gwas_df1, &gwas_df2, false)?;

    let output_path = format_rg_output_path(&gwas_path_1, &gwas_path_2, &output_root);
    let progress_path = output_path.with_extension("progress.txt");

    let result = solve_cors_wrapper(
        &input_data.tagging,
        &input_data.gwas_sumstats1,
        &input_data.gwas_sumstats2,
        &input_data.category_values,
        &tag_info.category_info.ssums,
        progress_path.to_str().unwrap(),
        None,
    );

    let partitions = format_genetic_correlation(&result);

    write_results(&output_path, &partitions)?;
    drop(permit);

    Ok(())
}

pub async fn rg_aligned(
    tagging: &[f64],
    category_values: &[Vec<f64>],
    category_contribs: &[Vec<f64>],
    gwas_path_1: Arc<PathBuf>,
    gwas_path_2: Arc<PathBuf>,
    output_root: Arc<PathBuf>,
    semaphore: Arc<Semaphore>,
) -> Result<()> {
    let permit = semaphore.acquire().await?;
    let path_1 = gwas_path_1.clone();
    let path_2 = gwas_path_2.clone();
    let gwas_1 = spawn_blocking(move || read_gwas_aligned(path_1.as_path()));
    let gwas_2 = spawn_blocking(move || read_gwas_aligned(path_2.as_path()));
    let gwas_1 = gwas_1.await??;
    let gwas_2 = gwas_2.await??;

    let output_path = format_rg_output_path(&gwas_path_1, &gwas_path_2, &output_root);
    let progress_path = output_path.with_extension("progress.txt");

    let result = solve_cors_wrapper(
        tagging,
        &gwas_1,
        &gwas_2,
        category_values,
        category_contribs,
        progress_path.to_str().unwrap(),
        None,
    );

    let partitions = format_genetic_correlation(&result);

    write_results(&output_path, &partitions)?;
    drop(permit);

    Ok(())
}
