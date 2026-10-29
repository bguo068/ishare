use super::{BcftkError, Result};
use error_stack::*;

use rust_htslib::bcf;
use std::path::{Path, PathBuf};

/// open a bcf/vcf reader from a path, or from stdin if no path is given
pub fn bcf_reader(input: Option<&PathBuf>) -> Result<bcf::Reader> {
    match input {
        None => {
            eprintln!("INFO: input argument not provided, reading from stdin ...");
            bcf::Reader::from_stdin()
                .change_context(BcftkError::Input)
                .attach("cannot create bcf reader from stdin")
        }
        Some(p) => bcf::Reader::from_path(p)
            .change_context(BcftkError::Input)
            .attach_with(|| format!("cannot create bcf reader from {}", p.display())),
    }
}

/// map bcftools-style output format string to (uncompressed, format)
pub fn parse_bcf_output_format(fmt: &str) -> Result<(bool, bcf::Format)> {
    let res = match fmt {
        "b" => (false, bcf::Format::Bcf),
        "u" => (true, bcf::Format::Bcf),
        "z" => (false, bcf::Format::Vcf),
        "v" => (true, bcf::Format::Vcf),
        _ => bail!(BcftkError::Input
            .into_report()
            .attach(format!("invalid bcf/vcf output format: {fmt}"))),
    };
    Ok(res)
}

/// create the parent folders of an output file if they do not exist
pub fn create_parent_dir(p: impl AsRef<Path>) -> Result<()> {
    if let Some(parent) = p.as_ref().parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .change_context(BcftkError::Output)
                .attach_with(|| format!("cannot create folder {}", parent.display()))?;
        }
    }
    Ok(())
}
