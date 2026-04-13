use ahash::AHashMap;
use itertools::Itertools;
use std::path::PathBuf;

use crate::args::Level;
use crate::args::SharingMetric;
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
    min_sharing: f64,
    min_denominator: f64,
    min_numerator: f64,
    aggregate_only: bool,
    chunk_size: usize,
    metric: SharingMetric,
    min_ac: u32,
    max_ac: u32,
}

struct ChunkResult {
    paires_vec: Vec<(u32, u32, f64)>, // id1, id2, sharing
    npairs: u32,
    running_sum: f64,
}

pub fn main_rvshare(args: &Commands) -> Result<()> {
    if let Commands::RvShare {
        rec,
        id,
        groups,
        level,
        sharing_type,
        min_sharing,
        min_denominator,
        min_numerator,
        output,
        aggregate_only,
        chunk_size,
        metric,
        min_ac,
        max_ac,
    } = args
    {
        let min_jaccard = min_sharing.unwrap_or(-1.0f64);
        let min_total = min_denominator.unwrap_or(0.0f64);
        let min_shared = min_numerator.unwrap_or(0.0f64);

        // use this struct to avoid passing them individually in function calls
        let cli = CliParameters {
            level: *level,
            sharing_type: *sharing_type,
            min_sharing: min_jaccard,
            min_numerator: min_shared,
            min_denominator: min_total,
            aggregate_only: *aggregate_only,
            chunk_size: *chunk_size,
            metric: *metric,
            min_ac: *min_ac,
            max_ac: *max_ac,
        };

        if groups.is_some() && !id.is_some() {
            eprintln!("WARN: when --groups is set, --id options are ignored");
        }

        eprintln!("read and concat rare genotype");
        let (mut records, inds) = read_and_concat_rare_genotypes(rec, cli.min_ac, cli.max_ac)?;

        let freq_map = if matches!(cli.metric, SharingMetric::GRM) {
            eprintln!("build allele frequency map");
            crate::utils::calc_allele_frequency(&mut records, inds.v().len() * 2)
                .change_context(GtencodeError::Input)?
        } else {
            AHashMap::new()
        };

        eprintln!("sort_records according sharing level");
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

        eprintln!("prep groups mapping from group name to a vector of Ids in that group");
        let group_map = prep_groups(id, groups, level, &inds)?;

        // process id pairs in each group pair
        let mut id_pairs = vec![];
        for (grp1, ids1) in group_map.iter() {
            for (grp2, ids2) in group_map.iter() {
                // avoid repeated calculation such as grp1-grp2 and then grp2- rp1
                if grp1 < grp2 {
                    continue;
                }
                eprintln!("process id pairs in each group pair- {grp1} and {grp2}");
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
                let res = process_group_pair(&records, id_pairs.as_slice(), &cli, &freq_map)?;

                // write result for a group pair
                write_results_for_a_group_pair(
                    (grp1, grp2),
                    (ids1, ids2),
                    res,
                    &mut aggfile,
                    cli.aggregate_only,
                    cli.level,
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
    freq_map: &AHashMap<(u32, u8), f64>,
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
                let mut denominator1 = 0.0f64;
                let mut denominator2 = 0.0f64;
                let mut numerator: f64 = 0.0f64;

                match cli.metric {
                    SharingMetric::Jaccard => match cli.level {
                        Level::IndividualLevel => {
                            for (_, n1, n2) in
                                records.iter_individual_pair_pos_allele_count(id1, id2)
                            {
                                numerator += n1.min(n2) as f64;
                                denominator1 += n1.max(n2) as f64;
                            }
                        }
                        Level::HaplotypeLevel => {
                            for (_, n1, n2) in records.iter_genome_pair_pos_allele_count(id1, id2) {
                                numerator += n1.min(n2) as f64;
                                denominator1 += n1.max(n2) as f64;
                            }
                        }
                    },
                    SharingMetric::Cosine => {
                        match cli.level {
                            Level::IndividualLevel => {
                                for (_, n1, n2) in
                                    records.iter_individual_pair_pos_allele_count(id1, id2)
                                {
                                    numerator += n1 as f64 * n2 as f64;
                                    denominator1 += n1 as f64 * n1 as f64;
                                    denominator2 += n1 as f64 * n1 as f64;
                                }
                            }
                            Level::HaplotypeLevel => {
                                for (_, n1, n2) in
                                    records.iter_genome_pair_pos_allele_count(id1, id2)
                                {
                                    numerator += n1 as f64 * n2 as f64;
                                    denominator1 += n1 as f64 * n1 as f64;
                                    denominator2 += n1 as f64 * n1 as f64;
                                }
                            }
                        }
                        denominator2 = denominator2.sqrt();
                        denominator1 = denominator1.sqrt();
                    }
                    SharingMetric::GRM => match cli.level {
                        Level::IndividualLevel => {
                            for (pos_allele, n1, n2) in
                                records.iter_individual_pair_pos_allele_count(id1, id2)
                            {
                                let freq = freq_map[&pos_allele];
                                numerator += (n1 as f64 - 2.0 * freq) * (n2 as f64 - 2.0 * freq);
                                denominator1 += 2.0 * freq * (1.0 - freq);
                            }
                        }
                        Level::HaplotypeLevel => {
                            for (pos_allele, n1, n2) in
                                records.iter_genome_pair_pos_allele_count(id1, id2)
                            {
                                let freq = freq_map[&pos_allele];
                                numerator += (n1 as f64 - 2.0 * freq) * (n2 as f64 - 2.0 * freq);
                                denominator1 += 2.0 * freq * (1.0 - freq);
                            }
                        }
                    },
                }

                let denominator = denominator1 + denominator2;
                let sharing = numerator / denominator;
                let out = (id1, id2, sharing);
                if (denominator1 < cli.min_denominator)
                    || (numerator < cli.min_numerator)
                    || (sharing < cli.min_sharing)
                {
                    continue;
                }

                running_sum += sharing;

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
    ids: (&[u32], &[u32]),
    res: Vec<ChunkResult>,
    aggrefile: &mut std::io::BufWriter<std::fs::File>,
    aggregate_only: bool,
    level: Level,
    output: &Option<PathBuf>,
) -> Result<()> {
    let (grp1, grp2) = groups;
    let (ids1, ids2) = ids;
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
        for (id1, id2, sharing) in chunk_res.paires_vec {
            if output.is_none() {
                match level {
                    Level::HaplotypeLevel => {
                        println!("g1={id1}, g2={id2}, sharing: {sharing:.6}",);
                    }
                    Level::IndividualLevel => {
                        println!("ind1={id1}, ind2={id2}, sharing: {sharing:.6}",);
                    }
                }
            }
            // update the matrix
            resmat.set_by_names(id1, id2, sharing);
            if identifical {
                resmat.set_by_names(id2, id1, sharing);
            }
        }
    }
    let grand_mean = grand_running_sum / (grand_npairs as f64);
    use std::io::Write;

    writeln!(aggrefile, "{grp1}\t{grp2}\t{grand_mean}")
        .change_context(GtencodeError::Output)
        .attach("fail to write aggregate into file")?;

    if output.is_none() {
        println!("group1={grp1}, group2={grp2}, grand_mean={grand_mean}");
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
