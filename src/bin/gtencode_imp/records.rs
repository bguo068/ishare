use super::super::Commands;

use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, utils::path::from_prefix};

pub fn main_records(args: &Commands) {
    if let Commands::Records {
        rec,
        genome,
        pos,
        samples,
        out,
    } = args
    {
        // read records to file
        use std::time::Instant;
        let start = Instant::now();
        println!("# Loading genotype records ...");
        let mut records = GenotypeRecords::from_parquet_file(rec);
        let duration = start.elapsed();
        println!("# Loading Time : {:?}", duration);

        let mut choosen_genome: Vec<u32> = vec![];

        match genome {
            Some(g) => {
                choosen_genome.push(*g);
            }
            _ => {
                let inds = Individuals::from_parquet_file(
                    &from_prefix(&rec, "ind").expect("cannot get inds path"),
                );

                match samples.as_ref() {
                    Some(sfname) => {
                        for sample_name in std::fs::read_to_string(sfname)
                            .expect("cant open sample list file")
                            .trim()
                            .split('\n')
                        {
                            match inds.m().get(sample_name) {
                                Some(idx) => {
                                    choosen_genome.push(*idx as u32 * 2);
                                    choosen_genome.push(*idx as u32 * 2 + 1);
                                }
                                None => {
                                    panic!("sample name invalid");
                                }
                            }
                        }
                    }
                    None => (0..(inds.m().len() as u32)).for_each(|i| {
                        choosen_genome.push(i * 2);
                        choosen_genome.push(i * 2 + 1)
                    }),
                }
            }
        };

        choosen_genome.sort();

        records
            .subset_by_genomes(&choosen_genome.as_slice())
            .unwrap();

        match pos {
            Some(p) => {
                records.records_mut().retain(|x| x.get_position() == *p);
            }
            _ => {}
        };

        match out {
            Some(out) => {
                println!("output records counts: {}", records.records().len());
                records.into_parquet_file(&from_prefix(out, "rec").unwrap());
            }
            None => records.records().iter().for_each(|r| {
                print!("{:?}\n", r);
            }),
        }
    }
}
