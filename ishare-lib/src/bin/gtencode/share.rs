use super::{GtencodeError, Result};
use error_stack::*;

use super::Commands;
use ishare::genotype::rare::GenotypeRecords;

pub fn main_share(args: &Commands) -> Result<()> {
    if let Commands::Share { rec, a, b } = args {
        let records =
            GenotypeRecords::from_parquet_file(rec).change_context(GtencodeError::Input)?;
        for (pos, gt1, gt2) in records.iter_genome_pair_genotypes(*a, *b) {
            println!("pos={pos}, allele_a={gt1:?}, allelle_b={gt2:?}");
        }
    }
    Ok(())
}
