use anyhow::{Context, Result};
use csv::ReaderBuilder;
use ndarray::prelude::*;
use polars::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

#[derive(Debug, Clone)]
pub struct CategoryInfo {
    pub n_variants: usize,
    pub n_categories: usize,
    pub names: Vec<String>,
    pub ssums: Vec<Vec<f64>>,
}

fn header_line_to_category_names(header_line: &str) -> Vec<String> {
    header_line
        .split_whitespace()
        .skip(9)
        .map(|x| x.to_string())
        .collect::<Vec<String>>()
}

fn final_lines_to_ssums(contrib_lines: &[String], final_line: &str) -> Vec<Vec<f64>> {
    // Contrib lines look like below. We want [1.2, 3.4, 5.6] etc.
    // "The relative contribution of the Base to each category 1.2 3.4 5.6" etc.
    let extract_contrib = |line: &str| {
        line.split_whitespace()
            .skip(9)
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<f64>>()
    };

    let mut ssums = contrib_lines
        .iter()
        .map(|x| extract_contrib(x))
        .collect::<Vec<Vec<f64>>>();

    // Final line looks like below. We want n_total = 100.
    // "There are 100 reference 99 regression 98 heritability predictors 97 96 95"
    let n_total = final_line
        .split_whitespace()
        .nth(2)
        .unwrap()
        .parse::<f64>()
        .unwrap();

    let final_vals = final_line
        .split_whitespace()
        .skip(9)
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<f64>>();

    // LDAK adds two values to each ssums vector: n_snps and n_snps / n_total
    final_vals
        .iter()
        .enumerate()
        .for_each(|(i, x)| ssums[i].push(*x));

    final_vals
        .iter()
        .enumerate()
        .for_each(|(i, x)| ssums[i].push(*x / n_total));

    ssums
}

fn read_category_info(filename: &str) -> Result<CategoryInfo> {
    let file = File::open(filename)?;
    let buf_reader = BufReader::new(file);

    let mut header_line = String::new();
    let mut n_variants = 0;
    let mut contrib_lines = Vec::new();
    let mut final_line = String::new();

    for (i, row) in buf_reader.lines().enumerate() {
        let line = row?;

        if i == 0 {
            header_line = line;
        } else if line.starts_with("The") {
            if n_variants == 0 {
                // Number of rows after header before tag info lines
                n_variants = i - 1;
            }

            if line.starts_with("There") {
                final_line = line;
            } else {
                contrib_lines.push(line);
            }
        }
    }

    let category_names = header_line_to_category_names(&header_line);
    let ssums = final_lines_to_ssums(&contrib_lines, &final_line);

    Ok(CategoryInfo {
        n_variants,
        n_categories: category_names.len(),
        names: category_names,
        ssums,
    })
}

#[derive(Clone)]
pub struct TagInfo {
    pub df: DataFrame,
    pub category_info: CategoryInfo,

    pub predictor_order: Vec<String>,
    pub tag_vec: Vec<f64>,
    pub cat_vec: Vec<Vec<f64>>,
}

impl TagInfo {
    pub fn from_dataframe(df: DataFrame, category_info: CategoryInfo) -> Result<Self> {
        let predictor_order = df
            .column("Predictor")?
            .utf8()?
            .into_iter()
            .map(|x| x.expect("Could not convert Predictor column").to_string())
            .collect::<Vec<String>>();

        let tag_vec = df
            .column("Tagging")?
            .f64()?
            .into_iter()
            .collect::<Option<Vec<f64>>>()
            .context("Tagging column contains null values!")?;

        let cat_vec = df
            .select(category_info.names.iter())?
            .to_ndarray::<Float64Type>(IndexOrder::Fortran)?
            .axis_iter(Axis(1))
            .map(|x| x.to_vec())
            .collect::<Vec<Vec<f64>>>();

        Ok(Self {
            df,
            category_info,
            predictor_order,
            tag_vec,
            cat_vec,
        })
    }

    pub fn update_from_dataframe(&mut self, df: DataFrame) -> Result<()> {
        let predictor_order = df
            .column("Predictor")?
            .utf8()?
            .into_iter()
            .map(|x| x.expect("Could not convert Predictor column").to_string())
            .collect::<Vec<String>>();

        let tag_vec = df
            .column("Tagging")?
            .f64()?
            .into_iter()
            .collect::<Option<Vec<f64>>>()
            .context("Tagging column contains null values!")?;

        let cat_vec = df
            .select(self.category_info.names.iter())?
            .to_ndarray::<Float64Type>(IndexOrder::Fortran)?
            .axis_iter(Axis(1))
            .map(|x| x.to_vec())
            .collect::<Vec<Vec<f64>>>();

        self.df = df;
        self.predictor_order = predictor_order;
        self.tag_vec = tag_vec;
        self.cat_vec = cat_vec;

        Ok(())
    }
}

pub fn read_tagfile(filename: &str) -> Result<TagInfo> {
    let category_info = read_category_info(filename)?;

    let mut csv_reader = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b' ')
        .from_path(filename)?;

    let header = csv_reader.headers()?.clone();
    let mut schema = Schema::new();

    header.iter().for_each(|x| {
        schema.with_column(
            x.into(),
            match x {
                "Predictor" => DataType::Utf8,
                _ => DataType::Float64,
            },
        );
    });

    let df = CsvReader::from_path(filename)?
        .has_header(true)
        .with_separator(b' ')
        .with_n_rows(Some(category_info.n_variants))
        .with_ignore_errors(true)
        .with_dtypes(Some(Arc::new(schema)))
        .finish()
        .with_context(|| format!("Failed to read tagfile {}", filename))?;

    TagInfo::from_dataframe(df, category_info)
}
