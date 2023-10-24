use super::super::Commands;
use super::utils::*;
use ishare::genotype::rare::GenotypeRecords;
use ishare::indiv::Individuals;
use ishare::io::IntoParquet;
use ishare::share::mat::NamedMatrix;
use rayon::prelude::*;
use std::path::PathBuf;

pub fn main_cosine(args: &Commands) {
    if let Commands::Cosine {
        rec,
        genomes,
        lists,
        min_cosine,
        min_magnitude,
        min_dot_prod,
        output,
    } = args
    {
        let min_cosine = min_cosine.unwrap_or(-0.001f64);
        let min_magnitude = min_magnitude.unwrap_or(-0.001f64);
        let min_dot_prod = min_dot_prod.unwrap_or(0i32);

        let records = GenotypeRecords::from_parquet_file(rec);
        assert!(records.is_sorted_by_genome());

        let ind_file = rec.with_extension("ind");
        let _inds = Individuals::from_parquet_file(&ind_file);

        let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);

        // run in parallel and collect row results
        let res: Vec<(u32, u32, i32, i32, i32)> = pairs
            .into_par_iter()
            .map(|(genome1, genome2)| {
                let mut sum_a2: i32 = 0;
                let mut sum_b2: i32 = 0;
                let mut sum_ab: i32 = 0;
                for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                    let (a, b) = match (a, b) {
                        (Some(a), Some(b)) if a == b => (1, 1),
                        // different rare variants, similarity decrease
                        (Some(_), Some(_)) => (-1, 1),
                        (Some(_), None) => (1, 0),
                        (None, Some(_)) => (0, 1),
                        (None, None) => (0, 0),
                    };
                    sum_a2 += a * a;
                    sum_b2 += b * b;
                    sum_ab += a * b;
                }

                let dot_prod = sum_ab;
                let magnitude = ((sum_a2 as f64) * (sum_b2 as f64)).sqrt();
                let cosine = dot_prod as f64 / magnitude;

                let mut out = Some((genome1, genome2, sum_a2, sum_b2, sum_ab));
                if (dot_prod < min_dot_prod) || (magnitude < min_magnitude) || (cosine < min_cosine)
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

        for (g1, g2, sum_a2, sum_b2, sum_ab) in res {
            let dot_prod = sum_ab;
            let magnitude = ((sum_a2 as f64) * (sum_b2 as f64)).sqrt();
            let cosine = dot_prod as f64 / magnitude;

            if output.is_none() {
                println!(
                        "g1={g1}, g2={g2}, sum_a2={sum_a2}, sum_b2={sum_b2}, sum_ab={sum_ab}, cosine={cosine:.6}",
                    );
            }

            // update the matrix
            resmat.set_by_names(g1, g2, cosine);
            if identifical {
                resmat.set_by_names(g2, g1, cosine);
            }
        }

        // write matrix to files
        match output.as_ref() {
            Some(output) => {
                let dir = output.parent().unwrap().to_str().unwrap();
                let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                if !filename.ends_with(".cos") {
                    filename.push_str(".cos");
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
