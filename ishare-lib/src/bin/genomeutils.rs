#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use std::{backtrace::Backtrace, path::PathBuf};

use clap::{Parser, Subcommand};
use ishare::{genome::Genome, utils::error::show_snafu_error};
use snafu::prelude::*;

#[derive(Debug, Snafu)]
enum Error {
    // #[snafu(transparent)]
    Genome {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genome::Error,
    },
    // local
    MissBothBppercmAndRate {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
}
type Result<T> = std::result::Result<T, Error>;
fn main() {
    if let Err(e) = main_entry() {
        show_snafu_error(e);
        std::process::exit(-1);
    }
}

#[derive(Parser, Debug)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Convert genome from text files to a binary file
    #[command(after_help=concat!(
        "Example: \n",
        "\tgenomeutils to-binary-file --from-toml tmp_genome.toml --to-bin tmp_genome.bin\n",
        "\n\tNOTE: an example of the genome.toml file can be generated with the ",
        "genomeutils generate-const-rr subcommand"
    ))]
    ToBinaryFile {
        /// Input genome TOML with the corrected gmaps paths
        #[arg(short = 't', long)]
        from_toml: PathBuf,

        /// Output binary file; if not set, it will be inferred from the --from-toml argument
        #[arg(short = 'b', long)]
        to_bin: Option<PathBuf>,
    },

    /// Convert genome from a binary file to text files
    ToTextFiles {
        /// Input binary file
        #[arg(short = 'b', long)]
        from_bin: PathBuf,

        /// Input genome TOML; if not set, it will be inferred from the --from-toml argument
        #[arg(short = 't', long)]
        to_toml: Option<PathBuf>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(short = 'x', long, default_value = "gmap/map")]
        to_map_prefix: PathBuf,
    },

    /// Generate genome from names of built-in genomes
    GenerateByName {
        /// Name of the built-in genome
        #[arg(short = 'n', long, required = true)]
        name: ishare::genome::BuiltinGenome,

        /// Output binary file
        #[arg(short = 'b', long, group = "generate_by_name_out_format")]
        to_bin: Option<PathBuf>,

        /// Output genome TOML
        #[arg(
            short = 't',
            long,
            group = "generate_by_name_out_format",
            group = "generate_by_name_to_toml"
        )]
        to_toml: Option<PathBuf>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(
            short = 'x',
            long,
            default_value = "gmap/map",
            requires = "generate_by_name_to_toml"
        )]
        to_map_prefix: PathBuf,
    },

    /// Generate genome using a constant recombination rate
    #[command(after_help = concat!(
        "Example:\n",
        "\tgenomeutils generate-const-rr \\\n",
        "\t\t--genome-name ex_genome  \\\n",
        "\t\t--chrom-size 100000000  \\\n",
        "\t\t--chrom-size 200000000  \\\n",
        "\t\t--chrom-name chr1  \\\n",
        "\t\t--chrom-name chr2  \\\n",
        "\t\t--bp-per-cm 100000000  \\\n",
        "\t\t--to-toml ex_genome.toml"))]
    GenerateConstRR {
        #[arg(short = 'N', long, default_value = "const_rate_genome")]
        genome_name: String,

        /// chromosome size in bp, if there are multiple chromosomes, the
        /// argument can be used multiple time. Note the number and orders
        /// of chromosome sizes and names should be consistent; at least one
        /// chromosome should be specified
        #[arg(short = 'l', long, required = true)]
        chrom_size: Vec<u32>,

        /// chromosome name, if there are multiple chromosomes, the
        /// argument can be used multiple time. Note the number and orders
        /// of chromosome sizes and names should be consistent; at least one
        /// chromosome should be specified
        #[arg(short = 'n', long, required = true)]
        chrom_name: Vec<String>,

        /// Recombination rate per bp per generation, e.g., --rate 1e-8
        #[arg(short = 'r', long, group = "grp_rate", required = true)]
        rate: Option<f32>,

        /// Recombination rate per bp per generation, e.g., --bp-per-cm 1000000
        #[arg(short = 'R', long, group = "grp_rate")]
        bp_per_cm: Option<u32>,

        /// Output binary file
        #[arg(short = 'b', long, group = "generate_by_name_out_format")]
        to_bin: Option<PathBuf>,

        /// Output genome TOML
        #[arg(short = 't', long, group = "generate_by_name_out_format")]
        to_toml: Option<PathBuf>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(
            short = 'x',
            long,
            default_value = "gmap/map",
            requires = "generate_by_name_to_toml"
        )]
        to_map_prefix: PathBuf,
    },
}

fn main_entry() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::ToBinaryFile { from_toml, to_bin } => {
            let p = to_bin.unwrap_or(from_toml.with_extension("bin"));
            Genome::load_from_text_file(&from_toml)
                .context(GenomeSnafu)?
                .save_to_bincode_file(&p)
                .context(GenomeSnafu)?;
        }
        Commands::ToTextFiles {
            from_bin,
            to_toml,
            to_map_prefix,
        } => {
            let p = to_toml.unwrap_or(from_bin.with_extension("bin"));
            let mut genome = Genome::load_from_bincode_file(&from_bin).context(GenomeSnafu)?;
            genome
                .save_to_text_files(&p.to_string_lossy(), &to_map_prefix.to_string_lossy())
                .context(GenomeSnafu)?;
        }
        Commands::GenerateByName {
            name,
            to_bin,
            to_toml,
            to_map_prefix,
        } => {
            let mut genome = Genome::new_from_name(name).context(GenomeSnafu)?;
            if let Some(to_bin) = to_bin {
                genome.save_to_bincode_file(&to_bin).context(GenomeSnafu)?;
            } else if let Some(to_toml) = to_toml {
                genome
                    .save_to_text_files(
                        &to_toml.to_string_lossy(),
                        &to_map_prefix.to_string_lossy(),
                    )
                    .context(GenomeSnafu)?;
            } else {
                eprintln!("Error: one of --to-bin and --to-toml has to be set");
                std::process::exit(-1);
            }
        }
        Commands::GenerateConstRR {
            genome_name,
            chrom_size,
            chrom_name,
            rate,
            bp_per_cm,
            to_bin,
            to_toml,
            to_map_prefix,
        } => {
            assert_eq!(chrom_size.len(), chrom_name.len());

            let rate = match rate {
                Some(rate) => rate,
                None => 0.01 / bp_per_cm.context(MissBothBppercmAndRateSnafu)? as f32,
            };
            let mut genome = Genome::new_from_constant_recombination_rate(
                &genome_name,
                &chrom_size,
                &chrom_name,
                rate,
            )
            .context(GenomeSnafu)?;
            if let Some(to_bin) = to_bin {
                genome.save_to_bincode_file(&to_bin).context(GenomeSnafu)?;
            } else if let Some(to_toml) = to_toml {
                genome
                    .save_to_text_files(
                        &to_toml.to_string_lossy(),
                        &to_map_prefix.to_string_lossy(),
                    )
                    .context(GenomeSnafu)?;
            } else {
                eprintln!("Error: one of --to-bin and --to-toml has to be set");
                std::process::exit(-1);
            }
        }
    }
    Ok(())
}
