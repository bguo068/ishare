use ishare::genotype::rare::GenotypeRecords;

use super::super::Commands;

use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    // non-local
    #[snafu(transparent)]
    GenotypeRare {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genotype::rare::Error,
    },
}
type Result<T> = std::result::Result<T, Error>;

pub fn main_share(args: &Commands) -> Result<()> {
    if let Commands::Share { rec, a, b } = args {
        let records = GenotypeRecords::from_parquet_file(rec)?;
        for (pos, gt1, gt2) in records.iter_genome_pair_genotypes(*a, *b) {
            println!("pos={pos}, allele_a={gt1:?}, allelle_b={gt2:?}");
        }
    }
    Ok(())
}
