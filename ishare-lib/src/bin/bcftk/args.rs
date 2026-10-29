use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name="bcftk", author, version, about, long_about=None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Clone)]
pub enum Commands {
    /// Filter Bcf/Vcf by per-sample and/or per-site genotype nonmissing rates
    Filtmiss {
        /// input bcf or vcf files
        #[arg(short, long)]
        input: PathBuf,

        /// filter string for setting nonmissing rates either for individuals (samples) - `i` or for variants (sites) - `v`.
        /// eg. `-f 0.3:v -f 0.4:i` will first filter sites and keep those with
        /// site non-missing rate >=0.3, then filter samples and keep those with
        /// sample non-missing rate >=0.4. Note there can be multiple rounds such as
        /// `-f 0.3:v -f 0.4:i -f 0.6:v -f 0.7:i`
        #[arg(short, long)]
        filter: Vec<String>,

        /// output prefix. Kept sites and samples are written to
        /// `{prefix}_sites.txt` and `{prefix}_samples.txt`
        #[arg(short, long, default_value = "out")]
        output_prefix: String,
    },

    /// Calculate the within-host diversity statistic Fws from FMT/AD of biallelic sites
    Calcfws {
        /// Optional input bcf/vcf file path. If not provided, read input from stdin
        input: Option<PathBuf>,

        /// Option output fws table. If not provided, will use `out_fws.tsv` as filename
        #[arg(short, long, value_name = "FILE")]
        output: Option<PathBuf>,
    },

    /// Drop non-snp (indel/mnp) alleles from each record, recoding FMT/GT and FMT/AD accordingly
    Dropindel {
        /// Optional input bcf/vcf file path. If not provided, read input from stdin
        input: Option<PathBuf>,

        /// Optional output bcf/vcf file path. If not provided, write to stdout
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,

        /// Output fmt: z: vcf.gz, v: vcf, b: bcf compressed, u: bcf uncompressed
        #[arg(short = 'O', long, default_value = "b")]
        output_format: String,

        /// Optional: maximum total frequency of non-snp alleles per record.
        /// If a site has a combined non-snp allele frequency higher than this threshold
        /// the whole record will be filtered out.
        #[arg(short = 'I', long)]
        max_total_nonsnp_afreq: Option<f64>,
    },
}
