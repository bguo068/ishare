use ishare::genotype::rare::GenotypeRecords;

use super::super::Commands;

pub fn main_share(args: &Commands) {
    if let Commands::Share { rec, a, b } = args {
        let records = GenotypeRecords::from_parquet_file(rec);
        for (pos, gt1, gt2) in records.iter_genome_pair_genotypes(*a, *b) {
            println!("pos={}, allele_a={:?}, allelle_b={:?}", pos, gt1, gt2);
        }
    }
}
