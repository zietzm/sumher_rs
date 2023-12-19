use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

mod ffi;
mod hsq;
mod io;
mod util;

use crate::util::{compute_hsq_parallel, compute_rg_parallel, format_plink_sumstats};

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
    Format {
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

fn validate_shared_args(args: &SharedArgs) {
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
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Command::Gcov {
            shared_args,
            h2_tagfile,
            rg_tagfile,
        } => {
            validate_shared_args(&shared_args);
            let result = compute_hsq_parallel(
                &h2_tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                shared_args.n_threads,
            );
            match result {
                Ok(_) => println!("Success on heritability!"),
                Err(e) => println!("Error: {}", e),
            }

            let result = compute_rg_parallel(
                &rg_tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                shared_args.n_threads,
            );
            match result {
                Ok(_) => println!("Success on genetic correlation!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::H2 {
            shared_args,
            tagfile,
        } => {
            validate_shared_args(&shared_args);
            let result = compute_hsq_parallel(
                &tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                shared_args.n_threads,
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
            validate_shared_args(&shared_args);
            let result = compute_rg_parallel(
                &tagfile,
                &shared_args.gwas_results,
                &shared_args.output_root,
                shared_args.n_threads,
            );
            match result {
                Ok(_) => println!("Success on genetic correlation!"),
                Err(e) => println!("Error: {}", e),
            }
        }
        Command::Format {
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

    // let path = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/plink_white_british.q_100001_0.glm.linear");
    // let output_path = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/plink_white_british.q_100001_0.glm.linear.summaries");
    // let result = format_plink_sumstats(path, output_path);
    // match result {
    //     Ok(_) => println!("Success!"),
    //     Err(e) => println!("Error: {}", e),
    // }

    // let output_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak/test/");
    // println!("Outputting to {}", output_root.to_str().unwrap());
    //
    // let tagroot = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example/ldak");
    // let tagpath = tagroot.join("ldak.thin.hapmap.gbr.tagging");
    // // let tagpath = tagroot.join("bld.ldak.hapmap.gbr.tagging");
    //
    // let gwas_root = Path::new("/Users/zietzm/Documents/projects/sumher_rs/example");
    // let gwas_paths = glob(gwas_root.join("*.glm.linear.summaries").to_str().unwrap())
    //     .unwrap()
    //     .map(|x| x.unwrap())
    //     .collect::<Vec<_>>();
    //
    // let gwas_paths = gwas_paths.iter().cycle().take(25).collect::<Vec<_>>();
    //
    // println!("Found {} GWAS file(s)", gwas_paths.len());
    //
    // // let result = compute_hsq_parallel(&tagpath, &gwas_paths, output_root, 500);
    // // match result {
    // //     Ok(_) => println!("Success on heritability!"),
    // //     Err(e) => println!("Error: {}", e),
    // // }
    //
    // let result = compute_rg_parallel(&tagpath, &gwas_paths, output_root, 100);
    // match result {
    //     Ok(_) => println!("Success on genetic correlation!"),
    //     Err(e) => println!("Error: {}", e),
    // }
    // println!("Done!");
}
