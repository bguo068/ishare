use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name="gtencode", author, version, about, long_about=None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Clone)]
pub enum Commands {
    /// Encode genotype data from VCF file to tables or matrix (-m)
    Encode {
        /// Path to VCF (input)
        vcf: PathBuf,
        /// Path to genome info toml file (input)
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to Genotype record table (output)
        #[arg(short, long)]
        out: PathBuf,
        /// Set chunks size for parallelized parsing, if not set, the whole genome is treated as a single chunk
        #[arg(short = 'c', long)]
        parallel_chunksize_bp: Option<usize>,
        /// If set, encode genotype as matrix (not efficient for rare variants);
        /// otherwise encode genotype as table (efficient for rare variants)
        #[arg(short, long, default_value_t = false)]
        matrix: bool,
        /// For table encoding only keep variants with maf no larger than the threshold;
        /// for matrix encoding only keep variants with maf no less than the threshold;
        #[arg(short = 'T', long, default_value_t = 0.001)]
        threshold_maf: f64,
    },

    /// View encoded table records
    Records {
        /// Path to Genotype record table (input)
        rec: PathBuf,
        /// optional genome id, if not specified, show all genomes
        #[arg(short, long)]
        genome: Option<u32>,
        /// optional genome-wide position, if not specified, show all positions
        #[arg(short, long)]
        pos: Option<u32>,
    },
    /// View encoded matrix contents
    Matrix {
        mat: PathBuf,
        /// optional genome id, if not specified, show all genomes
        #[arg(short, long)]
        genomes: Option<Vec<u32>>,
        /// optional genome-wide position, if not specified, show all positions
        #[arg(short, long)]
        positions: Option<Vec<u32>>,
    },
    /// View encoded sites
    Sites {
        /// Path to Genotype record table (input)
        sit: PathBuf,
        /// optional genome-wide pos, if not specified, show all positions
        #[arg(short, long)]
        pos: Option<u32>,
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
    },
    /// View individuals
    Samples {
        /// Path to Individal file (input)
        ind: PathBuf,
        /// optional sample name if not species show all samples.
        #[arg(short, long)]
        sample: Option<String>,
        /// optional sample idx if not species show all samples.
        #[arg(short = 'I', long)]
        idx_sample: Option<usize>,
        #[arg(short = 'i', long)]
        idx_genome: Option<usize>,
    },
    /// View genotype sharing between two genome
    Share {
        /// Path to Genotype record table (input)
        rec: PathBuf,
        /// genome 1
        a: u32,
        /// genome 2
        b: u32,
    },

    /// Calculate pairwise similarity via Jaccard index
    Jaccard {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of jaccard value of a genome pair to be printed
        #[arg(short = 'j', long)]
        min_jaccard: Option<f64>,
        /// optional low threshold of total number of sites with rare alleles for a genome pair
        #[arg(short = 't', long)]
        min_total: Option<u32>,
        /// optional low threshold of total number of sites where a genome pair shares a rare allele
        #[arg(short = 's', long)]
        min_shared: Option<u32>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },

    /// Calculate pairwise similarity via Cosine
    // See definition in https://en.wikipedia.org/wiki/Cosine_similarity
    Cosine {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of consine similarity for a genome pair
        #[arg(short = 'c', long)]
        min_cosine: Option<f64>,
        /// optional low threshold of magnitude (denominator) for a genome pair
        #[arg(short = 'm', long)]
        min_magnitude: Option<f64>,
        /// optional low threshold of dot product (numerator) for a genome pair
        #[arg(short = 'p', long)]
        min_dot_prod: Option<i32>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
    /// Calculate pairwise similarity via GRM (GCTA formula)
    // See defintion in Eqn3 of GCTA paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/
    Grm {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of grm relationship for a genome pair
        #[arg(short = 'r', long)]
        min_grm_related: Option<f64>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
    /// Run binary trait SKAT-O test
    // see details in Lee et al 2012 AJHG: https://www.cell.com/ajhg/fulltext/S0002-9297(12)00316-3
    Skato {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// windows size: the number of SNPs within the window
        #[arg(long, default_value_t = 200)]
        window: usize,
        /// step size: the number of SNPs between two consecutive windows's center
        #[arg(long, default_value_t = 200)]
        step: usize,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
    },
}
