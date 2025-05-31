use super::super::Commands;
use ishare::indiv::Individuals;

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    // non-local
    #[snafu(transparent)]
    Individuals { source: ishare::indiv::Error },
    // local
}

type Result<T> = std::result::Result<T, Error>;

pub fn main_samples(args: &Commands) -> Result<()> {
    if let Commands::Samples {
        ind,
        sample,
        idx_sample,
        idx_genome,
    } = args
    {
        let inds = Individuals::from_parquet_file(ind)?;

        let print_sample = |i: usize, s: &str| {
            println!(
                "idx={i}, genome1={}, genome2={}, name={s}",
                i * 2,
                i * 2 + 1
            );
        };

        #[derive(Copy, Clone)]
        enum Selection {
            None,
            One(usize),
            All,
        }

        let mut sel = match sample {
            Some(s) => Selection::One(inds.m()[s]),
            None => Selection::All,
        };
        sel = match (idx_sample, sel) {
            (Some(i), Selection::One(ii)) if *i == ii => sel,
            (Some(i), Selection::One(ii)) if *i != ii => Selection::None,
            (Some(_i), Selection::None) => Selection::None,
            (Some(i), Selection::All) => Selection::One(*i),
            _ => sel,
        };
        sel = match (idx_genome, sel) {
            (Some(i), Selection::One(ii)) if (*i / 2) == ii => sel,
            (Some(i), Selection::One(ii)) if (*i / 2) != ii => Selection::None,
            (Some(_i), Selection::None) => Selection::None,
            (Some(i), Selection::All) => Selection::One(*i / 2),
            _ => sel,
        };
        match sel {
            Selection::None => {}
            Selection::All => {
                for (i, s) in inds.v().iter().enumerate() {
                    print_sample(i, s);
                }
            }
            Selection::One(i) => {
                let s = &inds.v()[i];
                print_sample(i, s);
            }
        }
    }
    Ok(())
}
