use ahash::AHashMap;
use itertools::Itertools;
use slice_group_by::GroupByMut;
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
        from_processed_records,
        to_processed_records,
        target_group_pairs,
    } = args
    {
        let min_jaccard = min_sharing.unwrap_or(-1.0f64);
        let min_total = min_denominator.unwrap_or(0.0f64);
        let min_shared = min_numerator.unwrap_or(0.0f64);
        let min_ac = min_ac.as_ref().unwrap_or(&0);
        let max_ac = max_ac.as_ref().unwrap_or(&u32::MAX);

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

        let (mut records, inds, freq_map) = if *from_processed_records {
            // read processed records
            let records =
                GenotypeRecords::from_parquet_file(&rec[0]).change_context(GtencodeError::Input)?;
            let ind_file = rec[0].with_extension("ind");
            let inds = ishare::indiv::Individuals::from_parquet_file(&ind_file)
                .change_context(GtencodeError::Input)
                .attach("fail to read individual file")?;
            let freq_map = if matches!(cli.metric, SharingMetric::GRM) {
                let mut records_for_freq =
                    GenotypeRecords::from_parquet_file(rec[0].with_extension(".rec2"))
                        .change_context(GtencodeError::Input)?;
                crate::utils::calc_allele_frequency(&mut records_for_freq, inds.v().len() * 2)
                    .change_context(GtencodeError::Input)?
            } else {
                AHashMap::new()
            };

            (records, inds, freq_map)
        } else {
            eprintln!("read and concat rare genotype");
            let (mut records, inds) = read_and_concat_rare_genotypes(rec, cli.min_ac, cli.max_ac)?;
            let freq_map = if matches!(cli.metric, SharingMetric::GRM) {
                eprintln!("build allele frequency map");
                let freq_map =
                    crate::utils::calc_allele_frequency(&mut records, inds.v().len() * 2)
                        .change_context(GtencodeError::Input)?;
                if *to_processed_records {
                    records
                        .clone()
                        .into_parquet_file(rec[0].with_extension(".rec2"))
                        .change_context(GtencodeError::Output)
                        .attach("fail to write records file for allele frequency calculation")?;
                }
                freq_map
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
            if *to_processed_records {
                let p = output
                    .as_ref()
                    .unwrap_or(&PathBuf::from("merged.rec"))
                    .with_extension("rec");
                records
                    .clone()
                    .into_parquet_file(&p)
                    .change_context(GtencodeError::Output)
                    .attach("fail to write merged records")?;
                inds.clone()
                    .into_parquet_file(p.with_extension("ind"))
                    .change_context(GtencodeError::Output)
                    .attach("fail to write ind file for merged records")?;
            }
            (records, inds, freq_map)
        };

        if *to_processed_records {
            eprintln!("successfully processed genotype records");
            return Ok(());
        }

        // create a file to write aggregate per group pairs
        let agg_file_path = output
            .as_ref()
            .unwrap_or(&PathBuf::from("gtencode_jaccard"))
            .with_extension("aggfile");
        let mut aggfile = std::fs::File::create(agg_file_path)
            // .map(std::io::BufWriter::new)
            .change_context(GtencodeError::Output)
            .attach("fail to create output file to write aggregates")?;

        eprintln!("prep groups mapping from group name to a vector of Ids in that group");
        let group_map = prep_groups(id, groups, level, &inds)?;

        let target_group_pairs =
            prep_target_group_pairs(&group_map, target_group_pairs, cli.sharing_type)?;

        // ensure level and sort status are consistent and consolidate records for only target samples
        eprintln!("consolidating records only for targeted individuals or genomes");
        let mut target_ids: Vec<u32> = vec![];
        {
            let mut tgrp = ahash::AHashSet::new();
            target_group_pairs.iter().for_each(|(grp1, grp2)| {
                tgrp.insert(grp1);
                tgrp.insert(grp2);
            });
            for grp in tgrp.iter() {
                target_ids.extend(group_map[*grp].iter());
            }
            target_ids.sort();
            eprintln!("target_ids counts: {}", target_ids.len());
        }

        match &level {
            Level::IndividualLevel => {
                ensure!(
                records.is_sorted_by_individual_position_allele(),
                GtencodeError::Input.into_report().attach("requested individual level sharing analysis but genotype records is not sorted by individual/position/allele")
                );

                records
                    .records_mut()
                    .linear_group_by_key_mut(|s| s.get_individual())
                    .merge_join_by(target_ids.iter(), |a, b| a[0].get_individual() < **b)
                    .for_each(|e| {
                        if let itertools::Either::Left(to_exclude) = e {
                            to_exclude.iter_mut().for_each(|rec| rec.set_sentinel());
                        }
                    });
                records.records_mut().retain(|rec| !rec.is_sentinel());
                eprintln!("record counts: {}", records.records().len());
            }
            Level::HaplotypeLevel => {
                ensure!(
                records.is_sorted_by_genome_position_allele(),
                GtencodeError::Input.into_report().attach("requested haplotype level sharing analysis but genotype records is not sorted by individual/position/allele")
                );
                records
                    .records_mut()
                    .linear_group_by_key_mut(|s| s.get_genome())
                    .merge_join_by(target_ids.iter(), |a, b| a[0].get_genome() < **b)
                    .for_each(|e| {
                        if let itertools::Either::Left(to_exclude) = e {
                            to_exclude.iter_mut().for_each(|rec| rec.set_sentinel());
                        }
                    });
                records.records_mut().retain(|rec| !rec.is_sentinel());
            }
        }

        // process id pairs in each group pair
        let mut id_pairs = vec![];
        for (grp1, grp2) in target_group_pairs.iter() {
            let ids1 = &group_map[grp1];
            let ids2 = &group_map[grp2];
            // avoid repeated calculation such as grp1-grp2 and then grp2- rp1
            if grp1 < grp2 {
                continue;
            }
            eprintln!("process group pair- GROUP1={grp1} and GROUP2={grp2}");
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
                        // order them
                        .map(|(a, b)| if a < b { (b, a) } else { (a, b) })
                        .collect_vec(),
                );
            };
            id_pairs.par_sort_unstable();

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

                npairs += 1;
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
    aggrefile: &mut std::fs::File,
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

    writeln!(aggrefile, "{grp1}\t{grp2}\t{grand_npairs}\t{grand_mean}")
        .change_context(GtencodeError::Output)
        .attach("fail to write aggregate into file")?;

    if output.is_none() {
        println!("group1={grp1}, group2={grp2}, npairs={grand_npairs}, grand_mean={grand_mean}");
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

fn prep_target_group_pairs(
    group_map: &std::collections::HashMap<String, Vec<u32>>,
    target_group_pairs_file: &Option<PathBuf>,
    sharing_type: SharingType,
) -> Result<Vec<(String, String)>> {
    let mut target_group_pair = vec![];

    if let Some(p) = target_group_pairs_file.as_ref() {
        for line in std::fs::read_to_string(p)
            .change_context(GtencodeError::Input)
            .attach("fail to read target group pair file")?
            .trim()
            .split("\n")
        {
            let mut fields = line.trim().split("\t");
            let mut grp1 = fields
                .next()
                .ok_or(GtencodeError::Input)
                .attach("fail to parse column 1 in target group file")?
                .to_owned();
            if !group_map.contains_key(&grp1) {
                eprintln!("{grp1} is not a valid group name, related group paris are ignored");
                continue;
            }
            let mut grp2 = fields
                .next()
                .ok_or(GtencodeError::Input)
                .attach("fail to parse column 2 in target group file")?
                .to_owned();
            if !group_map.contains_key(&grp2) {
                eprintln!("{grp2} is not a valid group name, related group paris are ignored");
                continue;
            }
            if grp1 < grp2 {
                std::mem::swap(&mut grp1, &mut grp2);
            }
            target_group_pair.push((grp1, grp2));
        }
    } else {
        for (grp1, _ids1) in group_map.iter() {
            for (grp2, _ids2) in group_map.iter() {
                // avoid repeated calculation such as grp1-grp2 and then grp2- rp1
                if grp1 < grp2 {
                    continue;
                }
                if grp1 == grp2 && !matches!(sharing_type, SharingType::BetweenGroup) {
                    target_group_pair.push((grp1.to_owned(), grp2.to_owned()))
                }
                if grp1 != grp2 && !matches!(sharing_type, SharingType::WithinGroup) {
                    target_group_pair.push((grp1.to_owned(), grp2.to_owned()))
                };
            }
        }
    }
    target_group_pair.sort();

    Ok(target_group_pair)
}
