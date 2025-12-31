use std::backtrace::Backtrace;
use std::str::Utf8Error;

use super::super::Commands;
use ishare::genome::GenomeInfo;
use ishare::site::Sites;
use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    FromUtf8 {
        // leaf
        source: Utf8Error,
        backtrace: Box<Option<Backtrace>>,
    },
    // #[snafu(transparent)]
    Site {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::site::Error,
    },
    // #[snafu(transparent)]
    Genome {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genome::Error,
    },
}
type Result<T> = std::result::Result<T, Error>;

pub fn main_sites(args: &Commands) -> Result<()> {
    if let Commands::Sites {
        sit,
        pos,
        genome_info,
    } = args
    {
        let ginfo = GenomeInfo::from_toml_file(genome_info).context(GenomeSnafu)?;

        let sites = Sites::from_parquet_file(sit).context(SiteSnafu)?;
        for i in 0..sites.len() {
            let (p, alleles) = sites.get_site_by_idx(i);
            match pos {
                Some(pos) if *pos != p => {
                    continue;
                }
                _ => {}
            }
            let alleles = std::str::from_utf8(alleles).context(FromUtf8Snafu)?;
            let (chrid, chrname, chrpos) = ginfo.to_chr_pos(p);
            println!(
                "gwpos={p}, chrid={chrid}, chrname={chrname}, pos={chrpos} alleles = {alleles}"
            );
        }
    }
    Ok(())
}
