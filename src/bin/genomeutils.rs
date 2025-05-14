use ishare::genome::*;

use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Convert genome from text files to a binary file
    ToBinaryFile {
        /// Input genome TOML with the corrected gmaps paths
        #[arg(short = 't', long)]
        from_toml: String,

        /// Output binary file; if not set, it will be inferred from the --from-toml argument
        #[arg(short = 'b', long)]
        to_bin: Option<String>,
    },

    /// Convert genome from a binary file to text files
    ToTextFiles {
        /// Input binary file
        #[arg(short = 'b', long)]
        from_bin: String,

        /// Input genome TOML; if not set, it will be inferred from the --from-toml argument
        #[arg(short = 't', long)]
        to_toml: Option<String>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(short = 'x', long, default_value = "gmap/map")]
        to_map_prefix: String,
    },

    /// Generate genome from names of built-in genomes
    GenerateByName {
        /// Name of the built-in genome
        #[arg(short = 'n', long, required = true)]
        name: ishare::genome::BuiltinGenome,

        /// Output binary file
        #[arg(short = 'b', long, group = "generate_by_name_out_format")]
        to_bin: Option<String>,

        /// Output genome TOML
        #[arg(
            short = 't',
            long,
            group = "generate_by_name_out_format",
            group = "generate_by_name_to_toml"
        )]
        to_toml: Option<String>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(
            short = 'x',
            long,
            default_value = "gmap/map",
            requires = "generate_by_name_to_toml"
        )]
        to_map_prefix: String,
    },

    /// Generate genome using a constant recombination rate
    GenerateConstRR {
        #[arg(short = 'N', long, default_value = "const_rate_genome")]
        genome_name: String,

        #[arg(short = 'l', long, required = true)]
        chrom_size: Vec<u32>,

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
        to_bin: Option<String>,

        /// Output genome TOML
        #[arg(
            short = 't',
            long,
            group = "generate_by_name_out_format",
            group = "generate_by_name_to_toml"
        )]
        to_toml: Option<String>,

        /// Prefix for the output genome map files; defaults to "gmap/map"
        #[arg(
            short = 'x',
            long,
            default_value = "gmap/map",
            requires = "generate_by_name_to_toml"
        )]
        to_map_prefix: String,
    },
}
fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::ToBinaryFile { from_toml, to_bin } => {
            let p = to_bin.unwrap_or({
                std::path::Path::new(&from_toml)
                    .with_extension("bin")
                    .to_string_lossy()
                    .to_string()
            });
            Genome::load_from_text_file(&from_toml).save_to_bincode_file(&p);
        }
        Commands::ToTextFiles {
            from_bin,
            to_toml,
            to_map_prefix,
        } => {
            let p = to_toml.unwrap_or({
                std::path::Path::new(&from_bin)
                    .with_extension("bin")
                    .to_string_lossy()
                    .to_string()
            });
            let mut genome = Genome::load_from_bincode_file(&from_bin);
            genome.set_gmap_path_prefix(&to_map_prefix);
            genome.save_to_text_files(&p);
        }
        Commands::GenerateByName {
            name,
            to_bin,
            to_toml,
            to_map_prefix,
        } => {
            let mut genome = Genome::new_from_name(name);
            if let Some(to_bin) = to_bin {
                genome.save_to_bincode_file(&to_bin);
            } else if let Some(to_toml) = to_toml {
                genome.set_gmap_path_prefix(&to_map_prefix);
                genome.save_to_text_files(&to_toml);
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
                None => 0.01 / bp_per_cm.unwrap() as f32,
            };
            let mut genome = Genome::new_from_constant_recombination_rate(
                &genome_name,
                &chrom_size,
                &chrom_name,
                rate,
            );
            if let Some(to_bin) = to_bin {
                genome.save_to_bincode_file(&to_bin);
            } else if let Some(to_toml) = to_toml {
                genome.set_gmap_path_prefix(&to_map_prefix);
                genome.save_to_text_files(&to_toml);
            } else {
                eprintln!("Error: one of --to-bin and --to-toml has to be set");
                std::process::exit(-1);
            }
        }
    }
}
