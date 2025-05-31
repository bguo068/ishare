#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use args::{Cli, Commands};
use clap::Parser;
use gtencode_imp::*;
use ishare::utils::error::show_snafu_error;
use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    GtencodeSampleDiff { source: samplediff::Error },
    #[snafu(transparent)]
    GtencodeRvibd { source: rvibd::Error },
    #[snafu(transparent)]
    Export { source: export::Error },
    #[snafu(transparent)]
    GtencodeSites { source: sites::Error },
    #[snafu(transparent)]
    GtencodeCosine { source: cosine::Error },
    #[snafu(transparent)]
    GtencodeEncode { source: encode::Error },
    #[snafu(transparent)]
    GtencodeRecords { source: records::Error },
    #[snafu(transparent)]
    GtencodeMatrix { source: matrix::Error },
    #[snafu(transparent)]
    GtencodeSamples { source: samples::Error },
    #[snafu(transparent)]
    GtencodeShare { source: share::Error },
    #[snafu(transparent)]
    GtencodeJaccard { source: jaccard::Error },
    #[snafu(transparent)]
    GtencodeGrm { source: grm::Error },
    #[cfg(feature = "skato")]
    #[snafu(transparent)]
    GtencodeSkato { source: skato::Error },
}

type Result<T> = std::result::Result<T, Error>;

pub fn main() {
    if let Err(e) = main_entry() {
        show_snafu_error(e);
        std::process::exit(-1);
    }
}

mod gtencode_imp;
fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => encode::main_encode(args)?,
            args @ Commands::Records { .. } => records::main_records(args)?,
            args @ Commands::Matrix { .. } => matrix::main_matrix(args)?,
            args @ Commands::Sites { .. } => sites::main_sites(args)?,
            args @ Commands::Samples { .. } => samples::main_samples(args)?,
            args @ Commands::Share { .. } => share::main_share(args)?,
            args @ Commands::Jaccard { .. } => jaccard::main_jaccard(args)?,
            args @ Commands::SampleDiff { .. } => samplediff::main_samplediff(args)?,
            args @ Commands::Cosine { .. } => cosine::main_cosine(args)?,
            args @ Commands::Grm { .. } => grm::main_grm(args)?,
            #[cfg(feature = "skato")]
            args @ Commands::Skato { .. } => skato::main_skato(args)?,
            args @ Commands::RvIBD { .. } => rvibd::main_rvibd(args)?,
            args @ Commands::Export { .. } => export::main_export(args)?,
        },
        None => {
            println!("\nUse '-h  or [subcommand] -h' to show help message");
        }
    }
    Ok(())
}
