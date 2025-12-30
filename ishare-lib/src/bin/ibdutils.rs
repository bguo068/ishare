#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use clap::Parser;
pub mod ibdutils_impl;
use ibdutils_impl::{args::*, *};
use ishare::utils::error::show_snafu_error;
use snafu::prelude::*;

#[derive(Debug, Snafu)]
#[snafu(visibility)]
pub enum Error {
    #[snafu(transparent)]
    IbdutilsEncode {
        #[snafu(source(from(encode::Error, Box::new)))]
        source: Box<encode::Error>,
    },
    #[snafu(transparent)]
    IbdutilsCompare {
        #[snafu(source(from(compare::Error, Box::new)))]
        source: Box<compare::Error>,
    },
    #[snafu(transparent)]
    IbdutilsUnrelated {
        #[snafu(source(from(unrelated::Error, Box::new)))]
        source: Box<unrelated::Error>,
    },
    #[snafu(transparent)]
    IbdutilsCoverage {
        #[snafu(source(from(coverage::Error, Box::new)))]
        source: Box<coverage::Error>,
    },
    #[cfg(feature = "plotibd")]
    #[snafu(transparent)]
    IbdutilsPlotibd { source: plotibd::Error },
}
type Result<T> = std::result::Result<T, Error>;

pub fn main() {
    if let Err(e) = main_entry() {
        show_snafu_error(e);
        std::process::exit(-1);
    }
}

fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => encode::main_encode(args)?,
            args @ Commands::Compare { .. } => compare::main_compare(args)?,
            #[cfg(feature = "plotibd")]
            args @ Commands::PlotIBD { .. } => plotibd::main_plotibd(args)?,
            args @ Commands::GetUnrelated { .. } => unrelated::main_unrelated(args)?,
            args @ Commands::Coverage { .. } => coverage::main_coverage(args)?,
        },

        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
    Ok(())
}
