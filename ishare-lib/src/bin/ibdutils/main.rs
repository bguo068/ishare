#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use error_stack::*;
pub type Result<T> = std::result::Result<T, Report<IbdUtilsError>>;

#[derive(Debug, thiserror::Error)]
pub enum IbdUtilsError {
    #[error("input error")]
    Input,
    #[error("output error")]
    Output,
    #[error("library error")]
    Library,
    #[error("zero number of chromosome")]
    ZeroNumChromosome,
}

pub mod args;
pub mod compare;
pub mod coverage;
pub mod encode;
#[cfg(feature = "plotibd")]
pub mod plotibd;
pub mod unrelated;
pub mod utils;

use clap::Parser;
// pub mod ibdutils_impl;
use args::*;

pub fn main() -> Result<()> {
    main_entry()
}

fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => encode::main_encode(args).attach("ibdutil encode")?,
            args @ Commands::Compare { .. } => {
                compare::main_compare(args).attach("ibdutil compare")?
            }
            #[cfg(feature = "plotibd")]
            args @ Commands::PlotIBD { .. } => {
                plotibd::main_plotibd(args).attach("ibdutil plot-ibd")?
            }
            args @ Commands::GetUnrelated { .. } => {
                unrelated::main_unrelated(args).attach("ibdutil get-unrelated")?
            }
            args @ Commands::Coverage { .. } => {
                coverage::main_coverage(args).attach("ibdutil coverage")?
            }
        },

        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
    Ok(())
}
