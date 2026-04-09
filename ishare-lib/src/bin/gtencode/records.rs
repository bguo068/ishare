use super::{GtencodeError, Result};
use error_stack::*;

use super::Commands;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, utils::path::from_prefix};

pub fn main_records(args: &Commands) -> Result<()> {
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
        let mut records =
            GenotypeRecords::from_parquet_file(rec).change_context(GtencodeError::Input)?;
        let duration = start.elapsed();
        println!("# Loading Time : {duration:?}");

        let mut choosen_genome: Vec<u32> = vec![];

        match genome {
            Some(g) => {
                choosen_genome.push(*g);
            }
            _ => {
                let inds = Individuals::from_parquet_file(
                    from_prefix(rec, "ind").change_context(GtencodeError::Input)?,
                )
                .change_context(GtencodeError::Input)?;

                match samples.as_ref() {
                    Some(sfname) => {
                        for sample_name in std::fs::read_to_string(sfname)
                            .change_context(GtencodeError::Input)?
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

        records = records
            .subset_by_genomes(choosen_genome.as_slice())
            .change_context(GtencodeError::Input)?;

        if let Some(p) = pos {
            records.records_mut().retain(|x| x.get_position() == *p);
        };

        match out {
            Some(out) => {
                println!("output records counts: {}", records.records().len());
                records
                    .into_parquet_file(
                        from_prefix(out, "rec").change_context(GtencodeError::Output)?,
                    )
                    .change_context(GtencodeError::Output)?;
            }
            None => records.records().iter().for_each(|r| {
                println!("{r:?}");
            }),
        }
    }
    Ok(())
}
