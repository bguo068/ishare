#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

use ishare::{
    genome::GenomeInfo,
    gmap::GeneticMap,
    indiv::Individuals,
    rfmix::{asibd::ASIBDSet, fb::*},
    share::ibd::ibdset::IbdSet,
    utils::error::show_snafu_error,
};
use std::{
    fs::File,
    io::{read_to_string, BufWriter},
    path::{Path, PathBuf},
};

use clap::Parser;

#[derive(Parser)]
#[command(name="asibd", author, version, about, long_about = None)]
struct Cli {
    /// path to rfmix2 fb.tsv file
    #[arg(short = 'f', long, required = true)]
    rfmix2_fb: PathBuf,

    /// path to hapibd ibd.gz file
    #[arg(short = 'i', long, required = true)]
    hapibd_ibd: PathBuf,

    /// path to genome toml file
    #[arg(short = 'g', long, required = true)]
    genome: PathBuf,

    /// path to sample list file. Only the listed samples will be analyzed. Each sample
    /// in the list should be present in the rfmix2 fb.csv file header
    #[arg(short = 's', long, required = true)]
    samples: PathBuf,

    /// path to a list of snp psitions (1-based). These positions should be
    /// present both in the query and the reference vcf files that are used to call
    /// Local ancestry. The file has no header; each line has
    /// two fields (seperated by tab): chroname name and snp bp pos
    #[arg(short = 'p', long, required = true)]
    positions: PathBuf,

    /// minimal length in cM of the input IBD segments that will be analyzed to
    /// generate ancestry specific IBD segments. If not set,
    /// all IBD segments will be analyzed.
    #[arg(short = 'l', long)]
    min_ibd_seg_cm: Option<f32>,

    /// minimal probability in fb.tsv files that will allow ancestry assignment.
    /// If probabilities for all ancestry for a haplotype at a give site are all
    /// below this threshold, ancestry assignment will be set to unknown.
    #[arg(short = 'P', long, default_value_t = 0.9)]
    min_prob: f32,

    /// size of memory buffer for reading the large fb.tsv file in Mb.
    /// Large buffer size might reduce the number of system calls for io over network
    /// and accelerate the file reading step
    #[arg(short = 'B', long, default_value_t = 100)]
    buffer_size_mb: usize,

    /// path to output file
    #[arg(short = 'p', long, required = true)]
    out: PathBuf,
}

type Result<T> = std::result::Result<T, Error>;

use snafu::prelude::*;
use std::backtrace::Backtrace;

#[derive(Debug, Snafu)]
#[snafu(visibility)]
pub enum Error {
    #[snafu(transparent)]
    Io {
        source: ishare::io::Error,
    },
    #[snafu(transparent)]
    Sites {
        source: ishare::site::Error,
    },
    #[snafu(transparent)]
    Indiv {
        source: ishare::indiv::Error,
    },
    #[snafu(transparent)]
    Genome {
        source: ishare::genome::Error,
    },
    #[snafu(transparent)]
    GenotypeRare {
        source: ishare::genotype::rare::Error,
    },
    #[snafu(transparent)]
    Vcf {
        source: ishare::vcf::Error,
    },
    #[snafu(transparent)]
    Gmap {
        source: ishare::gmap::Error,
    },
    #[snafu(transparent)]
    Ibd {
        source: ishare::share::ibd::Error,
    },
    #[snafu(transparent)]
    Asibd {
        source: ishare::rfmix::asibd::Error,
    },

    #[snafu(transparent)]
    StdIo {
        source: std::io::Error,
    },
    #[snafu(transparent)]
    RfmixFb {
        source: ishare::rfmix::fb::Error,
    },

    RecordNotSorted {
        backtrace: Option<Backtrace>,
    },
    ReadChrname,
    ReadChrpos,
}
fn main() {
    if let Err(e) = main_entry() {
        show_snafu_error(e);
        std::process::exit(-1);
    }
}

fn main_entry() -> Result<()> {
    let cli = Cli::parse();

    eprintln!("reading gnome.toml");
    let ginfo = GenomeInfo::from_toml_file(&cli.genome)?;
    eprintln!("reading genetic map files");
    let gmap = GeneticMap::from_genome_info(&ginfo)?;
    eprintln!("reading sample list");
    let (indivs, opt) = Individuals::from_txt_file(&cli.samples)?;
    assert!(opt.is_none(), "should use single column samples list");

    eprintln!("reading ibd file");
    let mut ibd = IbdSet::new(&gmap, &ginfo, &indivs);
    ibd.read_hapibd_file(&cli.hapibd_ibd, cli.min_ibd_seg_cm)?;
    ibd.sort_by_haplotypes();

    eprintln!("reading LA position list");
    let mut pos = read_chr_positions(&cli.positions, &ginfo)?;
    pos.sort();
    eprintln!("reading LA fb.tsv file");
    let (ancestry, la_set) = {
        let fb = FbMatrix::from_fb_csv(
            &cli.rfmix2_fb,
            &pos,
            &ginfo,
            &indivs,
            cli.min_prob,
            cli.buffer_size_mb,
        )?;
        let ancestry = fb.get_ancestries().to_owned();
        let la_set = LASet::from_fbmat(&fb);
        (ancestry, la_set)
    };

    eprintln!("inferring as-ibd file");
    let mut asibd = ASIBDSet::new(&gmap, &ginfo, &indivs, &ancestry);
    asibd.get_asibd_from_ibdsets_and_laset(&ibd, &la_set)?;

    eprintln!("writing as-ibd file");
    let out = File::create(&cli.out).map(BufWriter::new)?;
    asibd.flush(out)?;
    Ok(())
}

fn read_chr_positions(p: impl AsRef<Path>, ginfo: &GenomeInfo) -> Result<Vec<u32>> {
    let buf = File::open(p.as_ref()).map(read_to_string)??;
    let mut v = vec![];
    buf.trim().split("\n").try_for_each(|x| {
        let mut splits = x.split("\t");
        let chrname = splits.next().context(ReadChrnameSnafu)?;
        let chr_pos = splits
            .next()
            .context(ReadChrnameSnafu)?
            .parse::<u32>()
            .ok()
            .context(ReadChrnameSnafu)?
            - 1; // 0-based position
        let chrid = ginfo.idx[chrname];
        v.push(ginfo.to_gw_pos(chrid, chr_pos));
        Ok::<(), Error>(())
    })?;
    Ok(v)
}
