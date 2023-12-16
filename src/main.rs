use glob::glob;
use ndarray::Axis;
use polars::prelude::*;
use std::error::Error;
use std::path::Path;

mod ffi;
mod io;
use crate::ffi::solve_sums_wrapper;
use crate::io::gwas::read_gwas_result;
use crate::io::tagging::{read_tagfile, TagInfo};

fn estimate_heritability(tag_info: &TagInfo, gwas_path: &str) -> Result<(), Box<dyn Error>> {
    // 1. Read GWAS file
    let gwas_df = read_gwas_result(gwas_path)?;

    println!("Tag schema: {:?}", tag_info.df.schema());
    println!("GWAS schema: {:?}", gwas_df.schema());

    // 2. Process data (merge with tag_info and extract to raw arrays)
    let full_df = gwas_df.join(
        &tag_info.df,
        ["Predictor"],
        ["Predictor"],
        JoinArgs::new(JoinType::Inner),
    )?;

    let n_variants = full_df.height() as i32;

    let tagging = full_df
        .column("Tagging")?
        .f64()?
        .into_iter()
        .map(|x| match x {
            Some(x) => x,
            None => panic!("Tagging column contains null values!"),
        })
        .collect::<Vec<f64>>();

    let category_values = full_df
        .select(tag_info.category_info.category_names.iter())?
        .to_ndarray::<Float64Type>(IndexOrder::Fortran)?
        .axis_iter(Axis(1))
        .map(|x| x.to_vec())
        .collect::<Vec<Vec<f64>>>();

    let category_contribs = &tag_info.category_info.ssums;

    let sample_sizes = full_df
        .column("n")?
        .i64()?
        .into_iter()
        .map(|x| match x {
            Some(x) => x as f64,
            None => panic!("n column contains null values!"),
        })
        .collect::<Vec<f64>>();

    let chisq = full_df
        .column("Z")?
        .f64()?
        .into_iter()
        .map(|x| match x {
            Some(x) => x.powi(2),
            None => panic!("Z column contains null values!"),
        })
        .collect::<Vec<f64>>();

    // 3. Pass to sumher using FFI
    let _result = solve_sums_wrapper(
        &tagging,
        &category_values,
        category_contribs,
        &sample_sizes,
        &chisq,
        n_variants,
        tag_info.category_info.n_categories as i32,
        "progress.txt\x00",
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    );
    Ok(())
}

fn main() {
    let tagroot = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak");
    let tagpath = tagroot.join("ldak.thin.hapmap.gbr.tagging");
    // let tagpath = tagroot.join("ldak.thin.hapmap.gbr.tagging_TEST");
    // let tagpath = tagroot.join("bld.ldak.hapmap.gbr.tagging");

    if !tagpath.exists() {
        println!("File {} does not exist!", tagpath.to_str().unwrap());
        return;
    }
    let tag_info = read_tagfile(tagpath.to_str().unwrap()).unwrap();

    let gwas_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example");
    let gwas_paths = glob(gwas_root.join("*.glm.linear.summaries").to_str().unwrap())
        .unwrap()
        .map(|x| x.unwrap().to_str().unwrap().to_string())
        .collect::<Vec<_>>();

    println!("Found {} GWAS file(s)", gwas_paths.len());
    let gwas_path = gwas_paths[0].as_str();

    let result = estimate_heritability(&tag_info, gwas_path);
    match result {
        Ok(_) => println!("Success!"),
        Err(e) => println!("Error: {}", e),
    }
}
