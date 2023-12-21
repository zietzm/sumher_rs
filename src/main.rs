use clap::{Args, Parser, Subcommand};
use hsq::compute_h2;
use std::path::PathBuf;

mod ffi;
mod hsq;
mod io;
mod util;

// use crate::hsq::compute_rg_parallel;
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
    },
}

fn validate_shared_args(args: &SharedArgs) -> RuntimeSetup {
    if args.gwas_results.is_empty() {
        panic!("No GWAS results files specified");
    }

    // Check whether the directory of output_root exists
    if !args.output_root.parent().unwrap().exists() {
        panic!(
            "Output root {} does not exist",
            args.output_root.parent().unwrap().to_str().unwrap()
        );
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.n_threads)
        .build_global()
        .unwrap();

    RuntimeSetup::new(args.n_threads)
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Command::Gcov {
            shared_args,
            h2_tagfile,
            rg_tagfile,
        } => {
            let rt = validate_shared_args(&shared_args);
            let result = compute_h2(
                &h2_tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                &rt,
            );
            // let result = compute_hsq_parallel(
            //     &h2_tagfile,
            //     &shared_args.gwas_results,
            //     &shared_args.output_root,
            //     &rt,
            // );
            match result {
                Ok(_) => println!("Success on heritability!"),
                Err(e) => println!("Error: {}", e),
            }

            // let result = compute_rg_parallel(
            //     &rg_tagfile,
            //     &shared_args.gwas_results,
            //     &shared_args.output_root,
            //     &rt,
            // );
            // match result {
            //     Ok(_) => println!("Success on genetic correlation!"),
            //     Err(e) => println!("Error: {}", e),
            // }
        }
        Command::H2 {
            shared_args,
            tagfile,
        } => {
            let rt = validate_shared_args(&shared_args);
            let result = compute_h2(
                &tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                &rt,
            );
            match result {
                Ok(_) => println!("Success on heritability!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::Rg {
            shared_args,
            tagfile,
        } => {
            let rt = validate_shared_args(&shared_args);
            // let result = compute_rg_parallel(
            //     &tagfile,
            //     &shared_args.gwas_results,
            //     &shared_args.output_root,
            //     &rt,
            // );
            // match result {
            //     Ok(_) => println!("Success on genetic correlation!"),
            //     Err(e) => println!("Error: {}", e),
            // }
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
