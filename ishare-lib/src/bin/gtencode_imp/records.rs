use std::backtrace::Backtrace;

use super::super::Commands;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, utils::path::from_prefix};

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    GenotypeRare {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genotype::rare::Error,
    },
    // #[snafu(transparent)]
    Individual {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::indiv::Error,
    },
    // #[snafu(transparent)]
    UtilsPath {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::utils::path::Error,
    },

    // local
    StdIo {
        // leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
}
type Result<T> = std::result::Result<T, Error>;

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
        let mut records = GenotypeRecords::from_parquet_file(rec).context(GenotypeRareSnafu)?;
        let duration = start.elapsed();
        println!("# Loading Time : {duration:?}");

        let mut choosen_genome: Vec<u32> = vec![];

        match genome {
            Some(g) => {
                choosen_genome.push(*g);
            }
            _ => {
                let inds = Individuals::from_parquet_file(
                    from_prefix(rec, "ind").context(UtilsPathSnafu)?,
                )
                .context(IndividualSnafu)?;

                match samples.as_ref() {
                    Some(sfname) => {
                        for sample_name in std::fs::read_to_string(sfname)
                            .context(StdIoSnafu)?
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
            .subset_by_genomes(choosen_genome.as_slice())
            .context(GenotypeRareSnafu)?;

        if let Some(p) = pos {
            records.records_mut().retain(|x| x.get_position() == *p);
        };

        match out {
            Some(out) => {
                println!("output records counts: {}", records.records().len());
                records
                    .into_parquet_file(from_prefix(out, "rec").context(UtilsPathSnafu)?)
                    .context(GenotypeRareSnafu)?;
            }
            None => records.records().iter().for_each(|r| {
                println!("{r:?}");
            }),
        }
    }
    Ok(())
}
