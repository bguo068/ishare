use clap::Parser;

pub mod ibdutils_impl;

use ibdutils_impl::{args::*, *};

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => encode::main_encode(args),
            args @ Commands::Compare { .. } => compare::main_compare(args),
            #[cfg(feature = "plotibd")]
            args @ Commands::PlotIBD { .. } => plotibd::main_plotibd(args),
            args @ Commands::GetUnrelated { .. } => unrelated::main_unrelated(args),
            args @ Commands::Coverage {.. } => coverage::main_coverage(args),
        },

        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
}
