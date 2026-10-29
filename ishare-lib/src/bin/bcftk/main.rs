#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]
use error_stack::*;
pub type Result<T> = std::result::Result<T, Report<BcftkError>>;

#[derive(Debug, thiserror::Error)]
pub enum BcftkError {
    #[error("input error")]
    Input,
    #[error("output error")]
    Output,
    #[error("library error")]
    Library,
}

pub mod args;
pub mod calcfws;
pub mod dropindel;
pub mod filtmiss;
pub mod utils;

use args::{Cli, Commands};
use clap::Parser;

pub fn main() -> Result<()> {
    main_entry()
}

fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Some(c) => match c {
            args @ Commands::Filtmiss { .. } => {
                filtmiss::main_filtmiss(args).attach("bcftk filtmiss")?
            }
            args @ Commands::Calcfws { .. } => {
                calcfws::main_calcfws(args).attach("bcftk calcfws")?
            }
            args @ Commands::Dropindel { .. } => {
                dropindel::main_dropindel(args).attach("bcftk dropindel")?
            }
        },
        None => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
    Ok(())
}
