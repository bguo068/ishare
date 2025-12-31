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
    // #[snafu(transparent)]
    IbdutilsEncode {
        // non leaf
        #[snafu(source(from(encode::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<encode::Error>,
    },
    // #[snafu(transparent)]
    IbdutilsCompare {
        // non leaf
        #[snafu(source(from(compare::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<compare::Error>,
    },
    // #[snafu(transparent)]
    IbdutilsUnrelated {
        // non leaf
        #[snafu(source(from(unrelated::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<unrelated::Error>,
    },
    // #[snafu(transparent)]
    IbdutilsCoverage {
        // non leaf
        #[snafu(source(from(coverage::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<coverage::Error>,
    },
    #[cfg(feature = "plotibd")]
    // #[snafu(transparent)]
    IbdutilsPlotibd {
        // non leaf
        #[snafu(backtrace)]
        #[snafu(source(from(plotibd::Error, Box::new)))]
        source: Box<plotibd::Error>,
    },
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
            args @ Commands::Encode { .. } => {
                encode::main_encode(args).context(IbdutilsEncodeSnafu)?
            }
            args @ Commands::Compare { .. } => {
                compare::main_compare(args).context(IbdutilsCompareSnafu)?
            }
            #[cfg(feature = "plotibd")]
            args @ Commands::PlotIBD { .. } => {
                plotibd::main_plotibd(args).context(IbdutilsPlotibdSnafu)?
            }
            args @ Commands::GetUnrelated { .. } => {
                unrelated::main_unrelated(args).context(IbdutilsUnrelatedSnafu)?
            }
            args @ Commands::Coverage { .. } => {
                coverage::main_coverage(args).context(IbdutilsCoverageSnafu)?
            }
        },

        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
    Ok(())
}
