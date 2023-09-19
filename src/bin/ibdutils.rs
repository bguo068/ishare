use clap::Parser;

pub mod ibdutils_impl;

use ibdutils_impl::{args::*, *};

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Compare { .. } => compare::main_compare(args),
            args @ Commands::PlotIBD { .. } => plotibd::main_plotibd(args),
            args @ Commands::GetUnrelated { .. } => unrelated::main_unrelated(args),
        },

        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
}
