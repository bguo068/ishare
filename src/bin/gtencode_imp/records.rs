use super::super::Commands;

use ishare::genotype::rare::GenotypeRecords;

pub fn main_records(args: &Commands) {
    if let Commands::Records { rec, genome, pos } = args {
        // read records to file
        use std::time::Instant;
        let start = Instant::now();
        println!("# Loading genotype records ...");
        let records = GenotypeRecords::from_parquet_file(rec);
        let duration = start.elapsed();
        println!("# Loading Time : {:?}", duration);
        for r in records.records() {
            let mut to_print = true;
            match pos {
                Some(p) if r.get_position() != *p => {
                    to_print = false;
                }
                _ => {}
            };
            match genome {
                Some(g) if r.get_genome() != *g => {
                    to_print = false;
                }
                _ => {}
            };
            if to_print {
                println!("{:?}", r)
            }
        }
    }
}
