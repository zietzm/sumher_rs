use glob::glob;
use std::path::Path;

use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

mod ffi;
mod hsq;
mod io;
use crate::hsq::estimate_heritability;
use crate::io::tagging::read_tagfile;

fn main() {
    let output_path = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak/");

    let tagroot = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak");
    // let tagpath = tagroot.join("ldak.thin.hapmap.gbr.tagging");
    let tagpath = tagroot.join("bld.ldak.hapmap.gbr.tagging");

    if !tagpath.exists() {
        println!("File {} does not exist!", tagpath.to_str().unwrap());
        return;
    }
    let tag_info = read_tagfile(tagpath.to_str().unwrap()).unwrap();

    let gwas_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example");
    let gwas_paths = glob(gwas_root.join("*.glm.linear.summaries").to_str().unwrap())
        .unwrap()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();

    println!("Found {} GWAS file(s)", gwas_paths.len());

    gwas_paths
        .iter()
        .cycle()
        .take(10)
        .collect::<Vec<_>>()
        // .iter()
        .into_par_iter()
        .progress_count(10_u64)
        .for_each(|x| {
            let output_path = output_path
                .join(x.file_stem().unwrap())
                .with_extension("hers");

            let result = estimate_heritability(
                &tag_info,
                x.to_str().unwrap(),
                output_path.to_str().unwrap(),
            );
            match result {
                Ok(_) => {}
                Err(e) => println!("Error: {}", e),
            }
        });

    println!("Done!");
}
