use crate::ffi::solve_sums_wrapper;
use crate::ffi::SolveSumsResult;
use crate::io::gwas::read_gwas_result;
use crate::io::tagging::TagInfo;

use ndarray::Axis;
use polars::prelude::*;

use std::error::Error;

use serde::Serialize;
use serde::Serializer;

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

fn write_heritability(
    gwas_path: &str,
    results: &[HeritabilityPartition],
) -> Result<(), Box<dyn Error>> {
    // Write using csv and serde serialization
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(gwas_path)?;

    for partition in results {
        writer.serialize(partition)?;
    }

    writer.flush()?;
    Ok(())
}

pub fn estimate_heritability(
    tag_info: &TagInfo,
    gwas_path: &str,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    // 1. Read GWAS file
    let gwas_df = read_gwas_result(gwas_path)?;

    // 2. Process data (merge with tag_info and extract to raw arrays)
    let full_df = gwas_df
        .join(
            &tag_info.df,
            ["Predictor"],
            ["Predictor"],
            JoinArgs::new(JoinType::Inner),
        )?
        .lazy()
        .with_column(col("Z").pow(2).alias("chisq"))
        .collect()?;

    let tagging = full_df
        .column("Tagging")?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or("Tagging column contains null values!")?;

    let sample_sizes = full_df
        .column("n")?
        .cast(&DataType::Float64)?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or("n column contains null values!")?;

    let chisq = full_df
        .column("chisq")?
        .f64()?
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or("Z column contains null values!")?;

    let category_values = full_df
        .select(tag_info.category_info.category_names.iter())?
        .to_ndarray::<Float64Type>(IndexOrder::Fortran)?
        .axis_iter(Axis(1))
        .map(|x| x.to_vec())
        .collect::<Vec<Vec<f64>>>();

    // 3. Pass to sumher using FFI
    let result = solve_sums_wrapper(
        &tagging,
        &chisq,
        &sample_sizes,
        &category_values,
        &tag_info.category_info.ssums,
        "progress.txt\x00",
        None,
    );

    // 4. Format result
    let partitions = format_heritability(&result, &tag_info.category_info.category_names);

    // 5. Write result
    write_heritability(output_path, &partitions)?;

    Ok(())
}
