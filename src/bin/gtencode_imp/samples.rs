use super::super::Commands;
use ishare::indiv::Individuals;

pub fn main_samples(args: &Commands) {
    if let Commands::Samples {
        ind,
        sample,
        idx_sample,
        idx_genome,
    } = args
    {
        let inds = Individuals::from_parquet_file(ind);

        let print_sample = |i: usize, s: &str| {
            println!(
                "idx={i}, genome1={}, genome2={}, name={s}",
                i * 2,
                i * 2 + 1
            );
        };

        #[derive(Copy, Clone)]
        enum Selection {
            NONE,
            ONE(usize),
            ALL,
        }

        let mut sel = match sample {
            Some(s) => Selection::ONE(inds.m()[s]),
            None => Selection::ALL,
        };
        sel = match (idx_sample, sel) {
            (Some(i), Selection::ONE(ii)) if *i == ii => sel,
            (Some(i), Selection::ONE(ii)) if *i != ii => Selection::NONE,
            (Some(_i), Selection::NONE) => Selection::NONE,
            (Some(i), Selection::ALL) => Selection::ONE(*i),
            _ => sel,
        };
        sel = match (idx_genome, sel) {
            (Some(i), Selection::ONE(ii)) if (*i / 2) == ii => sel,
            (Some(i), Selection::ONE(ii)) if (*i / 2) != ii => Selection::NONE,
            (Some(_i), Selection::NONE) => Selection::NONE,
            (Some(i), Selection::ALL) => Selection::ONE(*i / 2),
            _ => sel,
        };
        match sel {
            Selection::NONE => {}
            Selection::ALL => {
                for (i, s) in inds.v().iter().enumerate() {
                    print_sample(i, s);
                }
            }
            Selection::ONE(i) => {
                let s = &inds.v()[i];
                print_sample(i, s);
            }
        }
    }
}
