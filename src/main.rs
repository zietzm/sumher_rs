use clap::{Args, Parser, Subcommand};
use glob::glob;
use std::path::PathBuf;

mod ffi;
mod hsq;
mod io;
mod util;

// use crate::hsq::compute_rg_parallel;
use crate::hsq::{compute_h2, compute_rg};
use crate::util::{format_plink_sumstats, RuntimeSetup};

#[derive(Parser, Debug)]
#[command(name = "sumher_rs")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Args)]
#[command(args_conflicts_with_subcommands = true)]
#[command(flatten_help = true)]
struct SharedArgs {
    /// Paths to GWAS summary statistics files (in LDAK format)
    #[arg(num_args(1..), short, long, required = true)]
    gwas_results: Vec<PathBuf>,

    /// Stem of the output files
    #[arg(short, long, required = true)]
    output_root: PathBuf,

    /// Number of threads to use
    #[arg(short, long, default_value_t = num_cpus::get(), env = "NUM_THREADS")]
    n_threads: usize,

    /// Skip checking alignment (not recommended unless you know that all files are aligned)
    #[arg(long, default_value = "false")]
    skip_alignment_check: bool,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Format Plink summary statistics files for use with LDAK
    Fmt {
        /// Path to the GWAS summary statistics file
        #[arg(short, long)]
        gwas_results: PathBuf,

        /// Path to the output file
        #[arg(short, long)]
        output_path: PathBuf,
    },
    /// Compute full genetic covariance matrix (h2 and rg)
    Gcov {
        #[command(flatten)]
        shared_args: SharedArgs,

        /// Path to the LDAK tag file for heritability (e.g. bld.ldak.hapmap.gbr.tagging)
        #[arg(long)]
        h2_tagfile: PathBuf,

        /// Path to the LDAK tag file for genetic correlation (e.g. ldak.thin.hapmap.gbr.tagging)
        #[arg(long)]
        rg_tagfile: PathBuf,

        /// Number of phenotypes to load into memory at once
        #[arg(short, long, default_value_t = 100)]
        chunk_size: usize,
    },
    /// Compute just heritabilities
    H2 {
        #[command(flatten)]
        shared_args: SharedArgs,

        /// Path to the LDAK tag file (e.g. bld.ldak.hapmap.gbr.tagging)
        #[arg(short, long)]
        tagfile: PathBuf,
    },

    /// Compute just genetic correlations
    Rg {
        #[command(flatten)]
        shared_args: SharedArgs,

        /// Path to the LDAK tag file (e.g. ldak.thin.hapmap.gbr.tagging)
        #[arg(short, long)]
        tagfile: PathBuf,

        /// Number of phenotypes to load into memory at once
        #[arg(short, long, default_value_t = 100)]
        chunk_size: usize,
    },
}

fn validate_shared_args(args: &mut SharedArgs) -> RuntimeSetup {
    if args.gwas_results.is_empty() {
        panic!("No GWAS results files specified");
    }

    // Check whether the GWAS results files exist
    let mut new_gwas_results = Vec::new();
    for f in &args.gwas_results {
        if !f.exists() {
            if let Ok(globbed) = glob(f.to_str().unwrap()) {
                let globbed: Vec<_> = globbed.collect();

                if globbed.is_empty() {
                    panic!("GWAS results file {} does not exist", f.to_str().unwrap());
                }

                new_gwas_results.extend(globbed.into_iter().map(|x| x.unwrap()));
            }
        } else {
            new_gwas_results.push(f.clone());
        }
    }
    args.gwas_results = new_gwas_results;

    // Check whether the directory of output_root exists
    if !args.output_root.parent().unwrap().exists() {
        panic!(
            "Output root {} does not exist",
            args.output_root.parent().unwrap().to_str().unwrap()
        );
    }

    RuntimeSetup::new(args.n_threads)
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Command::Gcov {
            mut shared_args,
            h2_tagfile,
            rg_tagfile,
            chunk_size,
        } => {
            let rt = validate_shared_args(&mut shared_args);
            let result = compute_h2(
                &h2_tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                &rt,
                shared_args.skip_alignment_check,
            );
            match result {
                Ok(_) => println!("Success on heritability!"),
                Err(e) => println!("Error: {}", e),
            }

            let rg_result = compute_rg(
                &rg_tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                chunk_size,
                &rt,
                shared_args.skip_alignment_check,
            );
            match rg_result {
                Ok(_) => println!("Success on genetic correlation!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::H2 {
            mut shared_args,
            tagfile,
        } => {
            let rt = validate_shared_args(&mut shared_args);
            let result = compute_h2(
                &tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                &rt,
                shared_args.skip_alignment_check,
            );
            match result {
                Ok(_) => println!("Success on heritability!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::Rg {
            mut shared_args,
            tagfile,
            chunk_size,
        } => {
            let rt = validate_shared_args(&mut shared_args);
            let result = compute_rg(
                &tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                chunk_size,
                &rt,
                shared_args.skip_alignment_check,
            );
            match result {
                Ok(_) => println!("Success on genetic correlation!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::Fmt {
            gwas_results,
            output_path,
        } => {
            let result = format_plink_sumstats(&gwas_results, &output_path);
            match result {
                Ok(_) => println!("Success!"),
                Err(e) => println!("Error: {}", e),
            }
        }
    }
}
