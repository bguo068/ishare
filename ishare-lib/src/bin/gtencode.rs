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
    GtencodeSampleDiff {
        // non leaf
        #[snafu(source(from(samplediff::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<samplediff::Error>,
    },
    #[snafu(transparent)]
    GtencodeRvibd {
        // non leaf
        #[snafu(source(from(rvibd::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<rvibd::Error>,
    },
    #[snafu(transparent)]
    Export {
        // non leaf
        #[snafu(source(from(export::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<export::Error>,
    },
    #[snafu(transparent)]
    GtencodeStes {
        // non leaf
        #[snafu(source(from(sites::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<sites::Error>,
    },
    #[snafu(transparent)]
    GtencodeCosine {
        // non leaf
        #[snafu(source(from(cosine::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<cosine::Error>,
    },
    #[snafu(transparent)]
    GtencodeEncode {
        // non leaf
        #[snafu(source(from(encode::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<encode::Error>,
    },
    #[snafu(transparent)]
    GtencodeRecords {
        // non leaf
        #[snafu(source(from(records::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<records::Error>,
    },
    #[snafu(transparent)]
    GtencodeMatrix {
        // non leaf
        #[snafu(backtrace)]
        source: matrix::Error,
    },
    #[snafu(transparent)]
    GtencodeSamples {
        // non leaf
        #[snafu(backtrace)]
        source: samples::Error,
    },
    #[snafu(transparent)]
    GtencodeShare {
        // non leaf
        #[snafu(backtrace)]
        source: share::Error,
    },
    #[snafu(transparent)]
    GtencodeJaccard {
        // non leaf
        #[snafu(source(from(jaccard::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<jaccard::Error>,
    },
    #[snafu(transparent)]
    GtencodeGrm {
        // non leaf
        #[snafu(source(from(grm::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<grm::Error>,
    },
    #[cfg(feature = "skato")]
    #[snafu(transparent)]
    GtencodeSkato {
        // non leaf
        #[snafu(backtrace)]
        source: skato::Error,
    },
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
