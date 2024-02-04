use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::Result;

#[derive(Clone, Debug)]
pub struct TagLine {
    pub predictor: String,
    pub tagging: f64,
    pub annotations: Vec<f64>,
}

impl TagLine {
    pub fn new(predictor: String, tagging: f64, annotations: Vec<f64>) -> Self {
        Self {
            predictor,
            tagging,
            annotations,
        }
    }

    pub fn from_str(line: &str) -> Result<Self> {
        let parts = line.split_whitespace().collect::<Vec<_>>();
        let predictor = parts[0].to_string();
        let tagging = parts[4].parse()?;
        let annotations = parts[9..].iter().map(|x| x.parse().unwrap()).collect();
        Ok(Self::new(predictor, tagging, annotations))
    }
}

#[derive(Clone, Debug)]
pub struct TagInfo {
    pub taglines: Vec<TagLine>,
    pub annotation_names: Vec<String>,
    pub ssums: Vec<Vec<f64>>,

    pub predictors: Vec<String>,
    pub tagging: Vec<f64>,
    pub annotations: Vec<Vec<f64>>,
}

impl TagInfo {
    pub fn new(
        taglines: Vec<TagLine>,
        annotation_names: Vec<String>,
        ssums: Vec<Vec<f64>>,
    ) -> Self {
        let predictors = taglines.iter().map(|x| x.predictor.clone()).collect();
        let tagging = taglines.iter().map(|x| x.tagging).collect();

        let n_annotations = taglines[0].annotations.len();
        let mut annotations = vec![vec![0.0; taglines.len()]; n_annotations];

        for (i, tagline) in taglines.iter().enumerate() {
            for (j, annotation) in tagline.annotations.iter().enumerate() {
                annotations[j][i] = *annotation;
            }
        }

        Self {
            taglines,
            annotation_names,
            ssums,
            predictors,
            tagging,
            annotations,
        }
    }

    pub fn from_file<P>(file: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let file = File::open(file).unwrap();
        let reader = BufReader::new(file);

        let mut header_line = String::new();
        let mut tag_lines = Vec::new();
        let mut contrib_lines = Vec::new();
        let mut final_line = String::new();

        let mut reached_end = false;

        for (i, row) in reader.lines().enumerate() {
            let line = row?;

            if i == 0 {
                header_line = line;
            } else if line.starts_with("The") {
                reached_end = true;
                if line.starts_with("There") {
                    final_line = line;
                } else {
                    contrib_lines.push(line);
                }
            } else if !reached_end {
                tag_lines.push(TagLine::from_str(&line)?);
            }
        }

        let category_names = header_line_to_category_names(&header_line);
        let ssums = final_lines_to_ssums(&contrib_lines, &final_line);

        Ok(Self::new(tag_lines, category_names, ssums))
    }
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
