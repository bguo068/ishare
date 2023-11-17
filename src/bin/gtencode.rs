use args::{Cli, Commands};
use clap::Parser;
use gtencode_imp::*;

mod gtencode_imp;
fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Encode { .. } => encode::main_encode(args),
            args @ Commands::Records { .. } => records::main_records(args),
            args @ Commands::Matrix { .. } => matrix::main_matrix(args),
            args @ Commands::Sites { .. } => sites::main_sites(args),
            args @ Commands::Samples { .. } => samples::main_samples(args),
            args @ Commands::Share { .. } => share::main_share(args),
            args @ Commands::Jaccard { .. } => jaccard::main_jaccard(args),
            args @ Commands::SampleDiff { .. } => samplediff::main_samplediff(args),
            args @ Commands::Cosine { .. } => cosine::main_cosine(args),
            args @ Commands::Grm { .. } => grm::main_grm(args),
            #[cfg(feature = "skato")]
            args @ Commands::Skato { .. } => skato::main_skato(args),
            args @ Commands::RvIBD { .. } => rvibd::main_rvibd(args),
            args @ Commands::Export { .. } => export::main_export(args),
        },
        None => {
            println!("\nUse '-h  or [subcommand] -h' to show help message");
        }
    }
}
