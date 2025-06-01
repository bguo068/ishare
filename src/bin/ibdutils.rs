use clap::Parser;

pub mod ibdutils_impl;

use ibdutils_impl::{args::*, *};

use ishare::{gmap, utils::error::show_snafu_error};
use snafu::Snafu;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    GenomeError { source: ishare::genome::Error },
    #[snafu(transparent)]
    GmapError { source: gmap::Error },
}
pub(crate) type Result<T> = std::result::Result<T, Error>;

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
