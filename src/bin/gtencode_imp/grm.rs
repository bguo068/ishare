use super::super::Commands;
use super::utils::*;
use ishare::genotype::rare::GenotypeRecords;
use ishare::indiv::Individuals;
use ishare::io::IntoParquet;
use ishare::share::mat::NamedMatrix;
use ishare::site::Sites;
use rayon::prelude::*;
use std::path::PathBuf;

pub fn main_grm(args: &Commands) {
    if let Commands::Grm {
        rec,
        genomes,
        lists,
        min_grm_related,
        output,
    } = args
    {
        let min_grm_related = match min_grm_related {
            Some(x) => *x,
            None => -0.001f64,
        };

        let mut records = GenotypeRecords::from_parquet_file(rec);
        let ind_file = rec.with_extension("ind");
        let _inds = Individuals::from_parquet_file(&ind_file);
        let sites_file = rec.with_extension("sit");
        let _sites = Sites::from_parquet_file(&sites_file);
        let freq_map = calc_allele_frequency(&mut records, _inds.v().len() * 2, _sites.len());

        records.sort_by_genome();
        assert!(records.is_sorted_by_genome());

        let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);
        let base_sum = calc_base_relationship(&freq_map);

        // run in parallel and collect row results
        let res: Vec<(u32, u32, f64)> = pairs
            .into_par_iter()
            .map(|(genome1, genome2)| {
                let mut sum = base_sum;
                for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                    let (a, b) = match (a, b) {
                        (Some(a), Some(b)) if a == b => (0.0, 0.0),
                        // different rare variants, similarity decrease
                        (Some(_), Some(_)) => continue, // FIXME
                        (Some(_), None) => (0.0, 2.0),
                        (None, Some(_)) => (2.0, 0.0),
                        (None, None) => continue,
                    };
                    match freq_map.get(&_pos) {
                        Some(p) => {
                            sum -= (2.0 - 2.0 * p) * (2.0 - 2.0 * p) / 2.0 / p / (1.0 - p);
                            sum += (a - 2.0 * p) * (b - 2.0 * p) / 2.0 / p / (1.0 - p);
                        }
                        None => {}
                    }
                }

                let n = freq_map.len() as f64;
                let relationship = sum / n;

                let mut out = Some((genome1, genome2, relationship));
                if relationship < min_grm_related {
                    out = None
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

        for (g1, g2, relationship) in res {
            if output.is_none() {
                println!("g1={g1}, g2={g2}, relationship={relationship:.6}",);
            }

            // update the matrix
            resmat.set_by_names(g1, g2, relationship);
            if identifical {
                resmat.set_by_names(g2, g1, relationship);
            }
        }

        // write matrix to files
        match output.as_ref() {
            Some(output) => {
                let dir = output.parent().unwrap().to_str().unwrap();
                let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                if !filename.ends_with(".grm") {
                    filename.push_str(".grm");
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
