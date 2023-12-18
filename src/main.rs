use glob::glob;
use std::path::Path;

mod ffi;
mod hsq;
mod io;
mod util;
use crate::util::format_plink_sumstats;
use crate::util::{compute_hsq_parallel, compute_rg_parallel};

fn main() {
    // let path = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/plink_white_british.q_100001_0.glm.linear");
    // let output_path = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/plink_white_british.q_100001_0.glm.linear.summaries");
    // let result = format_plink_sumstats(path, output_path);
    // match result {
    //     Ok(_) => println!("Success!"),
    //     Err(e) => println!("Error: {}", e),
    // }

    let output_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak/");
    println!("Outputting to {}", output_root.to_str().unwrap());

    let tagroot = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak");
    let tagpath = tagroot.join("ldak.thin.hapmap.gbr.tagging");
    // let tagpath = tagroot.join("bld.ldak.hapmap.gbr.tagging");

    let gwas_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example");
    let gwas_paths = glob(gwas_root.join("*.glm.linear.summaries").to_str().unwrap())
        .unwrap()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();

    // Duplicate for testing
    let gwas_paths = gwas_paths
        .iter()
        .cycle()
        .take(100)
        // .map(|x| x.to_path_buf())
        .collect::<Vec<_>>();

    println!("Found {} GWAS file(s)", gwas_paths.len());

    let result = compute_hsq_parallel(&tagpath, &gwas_paths, output_root, 500);
    match result {
        Ok(_) => println!("Success on heritability!"),
        Err(e) => println!("Error: {}", e),
    }

    let result = compute_rg_parallel(&tagpath, &gwas_paths, output_root, 200);
    match result {
        Ok(_) => println!("Success on genetic correlation!"),
        Err(e) => println!("Error: {}", e),
    }
    println!("Done!");
}
