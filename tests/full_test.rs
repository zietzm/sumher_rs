use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::Command;

use assert_cmd::prelude::*;

#[test]
fn cli_h2_neur() {
    compute_h2("neur");
}

#[test]
fn cli_h2_height() {
    compute_h2("height");
}

fn compute_h2(phenotype: &str) {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("tests/data/");

    let output_dir = tempfile::tempdir().unwrap();
    let output = output_dir.path().join("output.db");

    let mut cmd = Command::cargo_bin("sumher_rs").unwrap();
    let assert = cmd
        .arg("h2")
        .arg("-g")
        .arg(d.join(phenotype.to_string() + ".txt").to_str().unwrap())
        .arg("-t")
        .arg(d.join("ldak.thin.genotyped.gbr.tagging").to_str().unwrap())
        .arg("-o")
        .arg(output.to_str().unwrap());

    assert.assert().success();

    let desired = load_ldak_h2_results(phenotype);
    let results = load_sumher_rs_h2_results(&output);

    assert_eq!(
        results.phenotype, desired.phenotype,
        "Phenotypes do not match"
    );
    assert_eq!(
        results.component, desired.component,
        "Components do not match"
    );
    assert_eq!(results.h2, desired.h2, "h2 values do not match");
    assert_eq!(results.sd, desired.sd, "Standard deviations do not match");
}

struct H2Result {
    phenotype: Vec<String>,
    component: Vec<String>,
    h2: Vec<f64>,
    sd: Vec<f64>,
}

fn load_sumher_rs_h2_results(db_path: &Path) -> H2Result {
    let conn = rusqlite::Connection::open(db_path).unwrap();
    let mut stmt = conn.prepare("SELECT * FROM h2").unwrap();

    let mut results = H2Result {
        phenotype: Vec::new(),
        component: Vec::new(),
        h2: Vec::new(),
        sd: Vec::new(),
    };

    stmt.query_map([], |row| {
        Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
    })
    .unwrap()
    .for_each(|x| {
        let (phenotype, component, h2, sd) = x.unwrap();
        results.phenotype.push(phenotype);
        // We made a conscious descision not to rename "Base" to "Her_Base" in
        // the sumher_rs code, so we need to do it here.
        if component == "Base" {
            results.component.push("Her_Base".to_string());
        } else {
            results.component.push(component);
        }
        results.h2.push(h2);
        results.sd.push(sd);
    });

    results
}

fn load_ldak_h2_results(phenotype: &str) -> H2Result {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/")
        .join(phenotype.to_owned() + ".hers");

    // Open the file at path
    let file = std::fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);

    let mut results = H2Result {
        phenotype: Vec::new(),
        component: Vec::new(),
        h2: Vec::new(),
        sd: Vec::new(),
    };

    reader.lines().skip(1).for_each(|line| {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split_whitespace().collect();
        results.phenotype.push(phenotype.to_string() + ".txt");
        results.component.push(fields[0].to_string());
        results.h2.push(fields[1].parse::<f64>().unwrap());
        results.sd.push(fields[2].parse::<f64>().unwrap());
    });

    results
}

// #[test]
// fn compute_rg() {
//     let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
//     d.push("tests/data/");
//
//     let output_dir = tempfile::tempdir().unwrap();
//     let output = output_dir.path().join("output.db");
//
//     let mut cmd = Command::cargo_bin("sumher_rs").unwrap();
//     let assert = cmd
//         .arg("rg")
//         .arg("-g")
//         .arg(d.join("neur.txt").to_str().unwrap())
//         .arg(d.join("height.txt").to_str().unwrap())
//         .arg("-t")
//         .arg(d.join("ldak.thin.genotyped.gbr.tagging").to_str().unwrap())
//         .arg("-o")
//         .arg(output.to_str().unwrap());
//
//     assert.assert().success();
//
//     let desired_results = load_ldak_rg_results();
//     let results = load_sumher_rs_rg_results(&output);
//
//     assert_eq!(
//         results.phenotype, desired_results.phenotype,
//         "Phenotypes do not match"
//     );
//     assert_eq!(
//         results.component, desired_results.component,
//         "Components do not match"
//     );
//     assert_eq!(results.h2, desired_results.h2, "h2 values do not match");
//     assert_eq!(
//         results.sd, desired_results.sd,
//         "Standard deviations do not match"
//     );
// }
//
// struct RgResult {
//     phenotype: Vec<String>,
//     component: Vec<String>,
//     value: Vec<f64>,
//     sd: Vec<f64>,
// }
//
// fn load_ldak_rg_results() -> RgResult {
//     let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
//         .join("tests/data/")
//         .join("neur.height.cors");
//
//     // Open the file at path
//     let file = std::fs::File::open(path).unwrap();
//     let reader = std::io::BufReader::new(file);
//
//     let mut results = RgResult {
//         phenotype: Vec::new(),
//         component: Vec::new(),
//         value: Vec::new(),
//         sd: Vec::new(),
//     };
//
//     reader.lines().skip(1).for_each(|line| {
//         let line = line.unwrap();
//         let fields: Vec<&str> = line.split_whitespace().collect();
//         results.phenotype.push(phenotype.to_string() + ".txt");
//         results.component.push(fields[0].to_string());
//         results.h2.push(fields[1].parse::<f64>().unwrap());
//         results.sd.push(fields[2].parse::<f64>().unwrap());
//     });
//
//     results
// }
