use itertools::Itertools;
use std::path::PathBuf;

use crate::args::Level;
use crate::args::SharingType;
use crate::utils::prep_groups;
use crate::utils::read_and_concat_rare_genotypes;

use super::{GtencodeError, Result};
use error_stack::*;

use super::Commands;
use ishare::genotype::rare::GenotypeRecords;
use ishare::io::IntoParquet;
use ishare::share::mat::NamedMatrix;
use rayon::prelude::*;

struct CliParameters {
    level: Level,
    sharing_type: SharingType,
    min_jaccard: f64,
    min_total: u32,
    min_shared: u32,
    aggregate_only: bool,
    chunk_size: usize,
}

struct ChunkResult {
    paires_vec: Vec<(u32, u32, f64)>, // id1, id2, sharing
    npairs: u32,
    running_sum: f64,
}

pub fn main_jaccard(args: &Commands) -> Result<()> {
    if let Commands::Jaccard {
        rec,
        id,
        groups,
        level,
        sharing_type,
        min_jaccard,
        min_total,
        min_shared,
        output,
        aggregate_only,
        chunk_size,
    } = args
    {
        let min_jaccard = min_jaccard.unwrap_or(-1.0f64);
        let min_total = min_total.unwrap_or(0u32);
        let min_shared = min_shared.unwrap_or(0u32);

        // use this struct to avoid passing them individually in function calls
        let cli = CliParameters {
            level: *level,
            sharing_type: *sharing_type,
            min_jaccard,
            min_shared,
            min_total,
            aggregate_only: *aggregate_only,
            chunk_size: *chunk_size,
        };

        if groups.is_some() && !id.is_some() {
            eprintln!("WARN: when --groups is set, --id options are ignored");
        }

        let (mut records, inds) = read_and_concat_rare_genotypes(rec)?;

        // sort_records according sharing level
        match &level {
            Level::IndividualLevel => {
                records
                    .sort_by_individual_position_allele()
                    .change_context(GtencodeError::Library)
                    .attach("fail to sort genotype records by individual/position/allele")?;
            }
            Level::HaplotypeLevel => {
                records
                    .sort_by_genome_position_allele()
                    .change_context(GtencodeError::Library)
                    .attach("fail to sort genotype records by genome position allele")?;
            }
        }

        // create a file to write aggregate per group pairs
        let agg_file_path = output
            .as_ref()
            .unwrap_or(&PathBuf::from("gtencode_jaccard"))
            .with_extension("aggfile");
        let mut aggfile = std::fs::File::create(agg_file_path)
            .map(std::io::BufWriter::new)
            .change_context(GtencodeError::Output)
            .attach("fail to create output file to write aggregates")?;

        // prep groups mapping from group name to a vector of Ids in that group
        let group_map = prep_groups(id, groups, level, &inds)?;

        // process id pairs in each group pair
        let mut id_pairs = vec![];
        for (grp1, ids1) in group_map.iter() {
            for (grp2, ids2) in group_map.iter() {
                id_pairs.clear();
                if grp1 == grp2 && !matches!(cli.sharing_type, SharingType::BetweenGroup) {
                    id_pairs.extend(
                        ids1.iter()
                            .copied()
                            .cartesian_product(ids2.iter().copied())
                            // filter to avoid duplicated calculation when a and b are in the same group
                            .filter(|(a, b)| a > b),
                    );
                }
                if grp1 != grp2 && !matches!(cli.sharing_type, SharingType::WithinGroup) {
                    id_pairs.extend(
                        ids1.iter()
                            .copied()
                            .cartesian_product(ids2.iter().copied())
                            .collect_vec(),
                    );
                };

                // process group pairs
                let res = process_group_pair(&records, id_pairs.as_slice(), &cli)?;

                // write result for a group pair
                write_results_for_a_group_pair(
                    (grp1, grp2),
                    ids1,
                    ids2,
                    res,
                    &mut aggfile,
                    cli.aggregate_only,
                    output,
                )?;
            }
        }
    }
    Ok(())
}

fn process_group_pair(
    records: &GenotypeRecords,
    pairs: &[(u32, u32)],
    cli: &CliParameters,
) -> Result<Vec<ChunkResult>> {
    // run in parallel and collect row results
    let res: Vec<_> = pairs
        .par_chunks(cli.chunk_size)
        .map(|chunk| {
            let mut local_vec = Vec::<(u32, u32, f64)>::with_capacity(chunk.len());
            let mut npairs = 0;
            let mut running_sum = 0.0f64;

            for &(id1, id2) in chunk {
                npairs += 1;
                let mut total: u32 = 0;
                let mut shared: u32 = 0;

                match cli.level {
                    Level::IndividualLevel => {
                        for (_pos, allele1_opt, allele2_opt) in
                            records.iter_individual_pair_genotypes(id1, id2)
                        {
                            match (allele1_opt, allele2_opt) {
                                (Some(_), Some(_)) => {
                                    shared += 1;
                                    total += 1
                                }
                                (None, None) => {}
                                _ => {
                                    total += 1;
                                }
                            }
                        }
                    }
                    Level::HaplotypeLevel => {
                        for (_pos, a, b) in records.iter_genome_pair_genotypes(id1, id2) {
                            match (a, b) {
                                (Some(a), Some(b)) if a == b => {
                                    shared += 1;
                                    total += 1
                                }
                                (None, None) => {}
                                (_, _) => total += 1,
                            }
                        }
                    }
                }

                let out = (id1, id2, shared as f64 / total as f64);
                if (total < cli.min_total)
                    || (shared < cli.min_shared)
                    || ((shared as f64) / (total as f64) < cli.min_jaccard)
                {
                    continue;
                }

                running_sum += (shared as f64) / (total as f64);

                if !cli.aggregate_only {
                    local_vec.push(out);
                }
            }
            ChunkResult {
                paires_vec: local_vec,
                npairs,
                running_sum,
            }
        })
        .collect::<Vec<_>>();
    Ok(res)
}

fn write_results_for_a_group_pair(
    groups: (&str, &str),
    ids1: &[u32],
    ids2: &[u32],
    res: Vec<ChunkResult>,
    aggrefile: &mut std::io::BufWriter<std::fs::File>,
    aggregate_only: bool,
    output: &Option<PathBuf>,
) -> Result<()> {
    let (grp1, grp2) = groups;
    // for identicial sets, elements correponsing to lower matrix is not
    // calculated use this as an indicator to fill the the low part of
    // the matrix when updating the full jaccard matrix
    let identifical = grp1 == grp2;

    let mut resmat = NamedMatrix::new(ids1.to_vec(), ids2.to_vec());

    let mut grand_npairs = 0;
    let mut grand_running_sum = 0.0f64;
    for chunk_res in res {
        grand_npairs += chunk_res.npairs;
        grand_running_sum += chunk_res.running_sum;
        for (g1, g2, sharing) in chunk_res.paires_vec {
            if output.is_none() {
                println!("g1={g1}, g2={g2}, sharing: {sharing:.6}",);
            }
            // update the matrix
            resmat.set_by_names(g1, g2, sharing);
            if identifical {
                resmat.set_by_names(g2, g1, sharing);
            }
        }
    }
    let grand_mean = grand_running_sum / (grand_npairs as f64);
    use std::io::Write;

    writeln!(aggrefile, "{grp1}\t{grp2}\t{grand_mean}")
        .change_context(GtencodeError::Output)
        .attach("fail to write aggregate into file")?;

    if output.is_none() {
        println!("grand_mean={grand_mean}");
    }

    if !aggregate_only {
        // write matrix to files
        if let Some(output) = output.as_ref() {
            let p = output.with_extension("jac.avg.txt");
            std::fs::write(p, format!("{grand_mean}"))
                .change_context(GtencodeError::Output)
                .attach("error write jav.avg.txt")?;
            let p = output.with_extension(format!("{grp1}_{grp2}.jac"));
            println!("WARN: output option is specified, results are not printed on the screen, check file {p:?}");
            // println!("\n writing...");
            resmat
                .into_parquet(&p)
                .change_context(GtencodeError::Output)?
        }
    }
    Ok(())
}
