use super::{GtencodeError, Result};
use error_stack::*;

use super::Commands;
use ishare::genome::GenomeInfo;
use ishare::site::Sites;

pub fn main_sites(args: &Commands) -> Result<()> {
    if let Commands::Sites {
        sit,
        pos,
        genome_info,
    } = args
    {
        let ginfo = GenomeInfo::from_toml_file(genome_info).change_context(GtencodeError::Input)?;

        let sites = Sites::from_parquet_file(sit).change_context(GtencodeError::Input)?;
        for i in 0..sites.len() {
            let (p, alleles) = sites.get_site_by_idx(i);
            match pos {
                Some(pos) if *pos != p => {
                    continue;
                }
                _ => {}
            }
            let alleles = std::str::from_utf8(alleles).change_context(GtencodeError::Input)?;
            let (chrid, chrname, chrpos) = ginfo.to_chr_pos(p);
            println!(
                "gwpos={p}, chrid={chrid}, chrname={chrname}, pos={chrpos} alleles = {alleles}"
            );
        }
    }
    Ok(())
}
