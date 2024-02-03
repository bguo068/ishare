use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name="ibdutils", author, version, about, long_about=None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Encode txt format IBD to binary version IBD with regions
    /// of extreme high IBD sharing or low SNP density region be excluded.
    Encode {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt: String,
        /// IBD directory
        #[arg(short = 'i', long, required = true)]
        ibd_dir: PathBuf,
        /// genome regions with IBD coverage larger
        /// than `outlier_upper` * std above 3% trimmed mean will be removed
        #[arg(long, default_value_t = 10.0)]
        outlier_upper: f64,
        /// position list to calcualte snp density. If None,
        // the position list will be generated from IBD segments
        #[arg(long)]
        position_lst: Option<PathBuf>,
        /// 1 cM windows with less than min_snp SNPs will be excluded
        #[arg(long, default_value_t = 15)]
        min_snp: u32,
        /// threshold: sample pairs with flattend IBD > genome size * threshold
        /// are treated as highly related.
        #[arg(long, default_value_t = 2.0)]
        min_cm: f32,
        /// Path to output file (prefix)
        #[arg(short = 'o', long, default_value = "ibd_encode")]
        out_prefix: PathBuf,
    },
    /// Encode genotype data from VCF file to tables or matrix (-m)
    Compare {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst1: PathBuf,
        /// Path to sample list file for ibd set2
        #[arg(short = 'S', long, required = true)]
        sample_lst2: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt1: String,
        /// fmt of ibd set2, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'F', long, required = true)]
        fmt2: String,
        /// IBD directory 1
        #[arg(short = 'i', long, required = true)]
        ibd1_dir: PathBuf,
        /// IBD directory 2
        #[arg(short = 'I', long, required = true)]
        ibd2_dir: PathBuf,
        /// by default no segments will be filterred out. By setting to a positive number
        /// any segment with length less than this positive number will be removed.
        #[arg(long, default_value_t = -0.1)]
        min_cm: f64,
        // length bin starts. example: "3,4,6,10,18"
        #[arg(long, default_value = "3,4,6,10,18")]
        length_bin_starts: String,
        // If specified, additional overlapping analysis will be done
        // for each window of the given size along the genome.
        #[arg(long)]
        window_size_bp: Option<u32>,
        /// Use (Do not Ignore) haplotype information (i.e. not flattening before overlapping analysis) if true
        #[arg(long, default_value_t = false)]
        use_hap_overlap: bool,
        /// Use (Do not Ignore) haplotype information for total ibd calculation. If true, calculate
        /// haplotype-level total ibd; if false, calculate merged sample-level total ibd.
        #[arg(long, default_value_t = false)]
        use_hap_totibd: bool,
        /// output details
        #[arg(short = 'D', long, default_value_t = false)]
        write_details: bool,
        /// Path to sample list file
        #[arg(short = 'o', long, default_value = "ibd_cmp_res.txt")]
        out: PathBuf,
    },
    #[cfg(feature = "plotibd")]
    PlotIBD {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst1: PathBuf,
        /// Path to sample list file for ibd set2
        #[arg(short = 'S', long, required = true)]
        sample_lst2: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt1: String,
        /// fmt of ibd set2, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'F', long, required = true)]
        fmt2: String,
        /// IBD directory 1
        #[arg(short = 'i', long, required = true)]
        ibd1_dir: PathBuf,
        /// IBD directory 2
        #[arg(short = 'I', long, required = true)]
        ibd2_dir: PathBuf,
        /// Path to output file (prefix)
        #[arg(short = 'o', long, default_value = "ibd_cmp_res")]
        out: PathBuf,

        #[command(flatten)]
        sample1: Sample1,
        #[command(flatten)]
        sample2: Sample2,

        /// Server at the given port
        #[arg(long)]
        port: Option<u16>,
    },

    /// Get a list of unrelated sample names.
    GetUnrelated {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt: String,
        /// IBD directory
        #[arg(short = 'i', long, required = true)]
        ibd_dir: PathBuf,
        /// threshold: sample pairs with flattend IBD > genome size * threshold
        /// are treated as highly related.
        #[arg(short = 'T', long, default_value_t = 0.5)]
        theshold: f64,
        /// Path to output file (prefix)
        #[arg(short = 'o', long, default_value = "unrelated.txt")]
        out: PathBuf,
    },
    /// Calculate IBD coverage for given sampling points
    Coverage {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set
        #[arg(short = 's', long, required = true)]
        sample_lst: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt: String,
        /// IBD directory
        #[arg(short = 'i', long, required = true)]
        ibd_dir: PathBuf,
        /// by default no segments will be filterred out. By setting to a
        /// positive number any segment with length less than this positive
        /// number will be removed.
        #[arg(long, default_value_t = -0.1)]
        min_cm: f64,
        /// position in CM of the most left sampling point
        #[arg(long, default_value_t = 0.01)]
        start_cm: f64,
        /// step size between adjacent sampling points
        #[arg(long, default_value_t = 0.01)]
        step_cm: f64,
        /// Use (Do not Ignore) haplotype information (i.e. not flattening
        /// before coverage analysis) if true
        #[arg(long, default_value_t = false)]
        prevent_flatten: bool,
        /// Path to output file in CSV format
        #[arg(short = 'o', long, default_value = "ibdcov_res.csv")]
        out: PathBuf,
    },
}

#[derive(Args)]
#[group(required = true, multiple = false)]
pub struct Sample1 {
    /// sample1 of the a pair by name
    #[arg(short = 'n', long)]
    pub ind_name1: Option<String>,
    /// sample1 of the a pair by index
    #[arg(short = 'x', long)]
    pub ind_ix1: Option<u32>,
}

#[derive(Args)]
#[group(required = true, multiple = false)]
pub struct Sample2 {
    /// sample2 of the a pair by name
    #[arg(short = 'N', long)]
    pub ind_name2: Option<String>,
    /// sample2 of the a pair by index
    #[arg(short = 'X', long)]
    pub ind_ix2: Option<u32>,
}
