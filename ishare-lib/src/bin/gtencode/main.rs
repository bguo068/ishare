#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]
use error_stack::*;
type Result<T> = std::result::Result<T, Report<GtencodeError>>;

#[derive(Debug, thiserror::Error)]
pub enum GtencodeError {
    #[error("input error")]
    Input,
    #[error("output error")]
    Output,
    #[error("library error")]
    Library,
}

pub mod args;
pub mod encode;
pub mod export;
pub mod matrix;
pub mod records;
pub mod rvibd;
pub mod rvshare;
pub mod samplediff;
pub mod samples;
pub mod sites;
#[cfg(feature = "skato")]
pub mod skato;
pub mod utils;

use args::{Cli, Commands};
use clap::Parser;

pub fn main() -> Result<()> {
    main_entry()
}

// mod gtencode_imp;
fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => {
                encode::main_encode(args).attach("gtencode encode")?
            }
            args @ Commands::Records { .. } => {
                records::main_records(args).attach("gtencode recode")?
            }
            args @ Commands::Matrix { .. } => {
                matrix::main_matrix(args).attach("gtencode matrix")?
            }
            args @ Commands::Sites { .. } => sites::main_sites(args).attach("gtencode sites")?,
            args @ Commands::Samples { .. } => {
                samples::main_samples(args).attach("gtencode samples")?
            }
            args @ Commands::RvShare { .. } => {
                rvshare::main_rvshare(args).attach("gtencode rv-share")?
            }
            args @ Commands::SampleDiff { .. } => {
                samplediff::main_samplediff(args).attach("gtencode sample-diff")?
            }
            #[cfg(feature = "skato")]
            args @ Commands::Skato { .. } => skato::main_skato(args).attach("gtencode skato")?,
            args @ Commands::RvIBD { .. } => rvibd::main_rvibd(args).attach("gtencode rv-ibd")?,
            args @ Commands::Export { .. } => {
                export::main_export(args).attach("gtencode export")?
            }
        },
        None => {
            println!("\nUse '-h  or [subcommand] -h' to show help message");
        }
    }
    Ok(())
}
