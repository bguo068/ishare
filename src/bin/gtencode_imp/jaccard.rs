use super::super::Commands;
use super::utils::*;
use ishare::genotype::rare::GenotypeRecords;
use ishare::indiv::Individuals;
use ishare::io::IntoParquet;
use ishare::share::mat::NamedMatrix;
use rayon::prelude::*;
use std::path::PathBuf;

pub fn main_jaccard(args: &Commands) {
    if let Commands::Jaccard {
        rec,
        genomes,
        lists,
        min_jaccard,
        min_total,
        min_shared,
        output,
    } = args
    {
        let min_jaccard = match min_jaccard {
            Some(x) => *x,
            None => -1.0f64,
        };
        let min_total = match min_total {
            Some(x) => *x,
            None => 0u32,
        };
        let min_shared = match min_shared {
            Some(x) => *x,
            None => 0u32,
        };

        let records = GenotypeRecords::from_parquet_file(rec);
        assert!(records.is_sorted_by_genome());

        let ind_file = rec.with_extension("ind");
        let _inds = Individuals::from_parquet_file(&ind_file);

        let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);

        // run in parallel and collect row results
        let res: Vec<(u32, u32, u32, u32)> = pairs
            .into_par_iter()
            .map(|(genome1, genome2)| {
                let mut total: u32 = 0;
                let mut shared: u32 = 0;
                for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                    match (a, b) {
                        (Some(a), Some(b)) if a == b => {
                            shared += 1;
                            total += 1
                        }
                        (None, None) => {}
                        (_, _) => total += 1,
                    }
                }

                let mut out = Some((genome1, genome2, total, shared));
                if (total < min_total)
                    || (shared < min_shared)
                    || ((shared as f64) / (total as f64) < min_jaccard)
                {
                    out = None;
                }

                out
            })
            .flatten()
            .collect();

        // for identicial sets, elements correponsing to lower matrix is not
        // calculated use this as an indicator to fill the the low part of
        // the matrix when updating the full jaccard matrix
        let identifical = row_genomes == col_genomes;

        let mut resmat = NamedMatrix::new(row_genomes, col_genomes);

        for (g1, g2, total, shared) in res {
            let jaccard = ((shared as f32) / (total as f32) * 1e6) as u32;

            if output.is_none() {
                println!(
                    "g1={g1}, g2={g2}, total={total}, shared={shared}, jaccard*10^6={}",
                    jaccard
                );
            }

            // update the matrix
            resmat.set_by_names(g1, g2, jaccard);
            if identifical {
                resmat.set_by_names(g2, g1, jaccard);
            }
        }

        // write matrix to files
        match output.as_ref() {
            Some(output) => {
                let dir = output.parent().unwrap().to_str().unwrap();
                let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                if !filename.ends_with(".jac") {
                    filename.push_str(".jac");
                }
                let mut p = PathBuf::from(dir);
                p.push(filename);

                // let resmat0 = resmat.clone();
                println!("WARN: output option is specified, results are not printed on the screen, check file {:?}", p);
                // println!("\n writing...");
                resmat.into_parquet(&p)
            }
            None => {}
        }
    }
}
