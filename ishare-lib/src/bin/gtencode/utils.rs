use crate::args::Level;

use super::{GtencodeError, Result};
use error_stack::*;

use ahash::AHashMap;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals};
use itertools::Itertools;
use log::warn;
use slice_group_by::{GroupBy, GroupByMut};
use std::{
    collections::{HashMap, HashSet},
    path::{Path, PathBuf},
};

pub fn file_to_u32_vec(p: impl AsRef<Path>) -> Result<Vec<u32>> {
    let mut s = String::new();
    let mut reader = std::fs::File::open(p)
        .map(std::io::BufReader::new)
        .change_context(GtencodeError::Input)?;
    s.clear();
    use std::io::Read;
    reader
        .read_to_string(&mut s)
        .change_context(GtencodeError::Input)?;
    let mut v: Vec<u32> = s
        .trim()
        .split("\n")
        .map(str::parse)
        .collect::<std::result::Result<_, _>>()
        .change_context(GtencodeError::Input)?;
    assert!(!v.is_empty());
    v.sort();
    Ok(v)
}

pub fn prep_groups(
    ids: &Option<Vec<u32>>,
    group_file: &Option<PathBuf>,
    level: &crate::args::Level,
    inds: &Individuals,
) -> Result<HashMap<String, Vec<u32>>> {
    let group_map = if let Some(ids) = ids {
        // get unique ids
        let mut v = ids.to_vec();
        v.sort_unstable();
        v.dedup();
        if v.len() < 2 {
            bail!(GtencodeError::Input
                .into_report()
                .attach("--id is provide but the number of unique id is less than two"));
        }
        if group_file.is_some() {
            warn!("when --id is provided, --groups is ignored")
        }
        let mut group_map = HashMap::with_capacity(1);
        group_map.insert("idx".to_owned(), v);
        group_map
    } else {
        // use --groups options
        let group_file = group_file
            .as_ref()
            .ok_or(GtencodeError::Input)
            .attach("neither --id nor --group are provided")?;
        read_groups_file(group_file, inds, *level).attach("fail to read group file")?
    };
    Ok(group_map)
}

pub type PairInfo = (
    Vec<(u32, u32)>, // pairs
    Vec<u32>,        // row genomes
    Vec<u32>,        // col genomes
);

pub fn prep_pairs(
    records: &GenotypeRecords,
    genomes: &Option<Vec<u32>>,
    lists: &[PathBuf],
) -> Result<PairInfo> {
    // genomes vec
    let genomes: Vec<u32> = match genomes {
        Some(v) => v.clone(),
        None => Vec::new(),
    };

    // read genome lists from files

    let mut lst1 = Vec::<u32>::new();
    let mut lst2 = Vec::<u32>::new();

    match lists.len() {
        0 => {}
        1 => {
            // within subset sharing
            lst1.extend(file_to_u32_vec(&lists[0])?);
            lst2.extend_from_slice(lst1.as_slice());
        }
        2 => {
            // inter-subset sharing
            lst1.extend(file_to_u32_vec(&lists[0])?);
            lst2.extend(file_to_u32_vec(&lists[1])?);

            // ensure non-overlapping
            let mut set = HashSet::<u32>::new();
            set.extend(&lst1);
            let overlapping = lst2.iter().any(|x| set.contains(x));
            if overlapping {
                eprintln!("list1 and list2 are overlapping");
                std::process::exit(-1);
            }
        }
        _ => {
            eprintln!("--lists options can specied no more than 2 times");
            std::process::exit(-1);
        }
    }

    let min_gid = records
        .records()
        .first()
        .ok_or(GtencodeError::Library)
        .attach("Record is empty")?
        .get_genome();
    let max_gid = records
        .records()
        .last()
        .ok_or(GtencodeError::Library)
        .attach("Record is empty")?
        .get_genome();
    let mut row_genomes = Vec::<u32>::new();
    let mut col_genomes = Vec::<u32>::new();
    let pairs: Vec<(u32, u32)> = match (!genomes.is_empty(), !lst1.is_empty()) {
        (true, false) => {
            row_genomes.extend(genomes.iter());
            col_genomes.extend(genomes.iter());
            genomes
                .iter()
                .copied()
                .cartesian_product(genomes.iter().copied())
                .filter(|(a, b)| *a > *b)
                .collect()
        }
        (false, true) => {
            row_genomes.extend(&lst1);
            col_genomes.extend(&lst2);

            if lists.len() == 2 {
                // two non-overallping list. no need to filter
                lst1.iter()
                    .copied()
                    .cartesian_product(lst2.iter().copied())
                    .collect()
            } else {
                // two identical list need to filtering
                lst1.iter()
                    .copied()
                    .cartesian_product(lst2.iter().copied())
                    .filter(|(a, b)| *a > *b)
                    .collect()
            }
        }
        (true, true) => {
            panic!("--gnomes and --lists should not be specified at the same time");
        }
        (false, false) => {
            let n = (max_gid - min_gid) as usize;
            if n * n * 8 / 1024 / 1024 / 1024 > 30 {
                eprintln!("too many pairs! consider using --lists option to restrict number of pairs for each run!");
                std::process::exit(-1);
            }
            row_genomes.extend(min_gid..max_gid);
            col_genomes.extend(min_gid..max_gid);
            (min_gid..max_gid)
                .cartesian_product(min_gid..max_gid)
                .filter(|(a, b)| *a > *b)
                .collect()
        }
    };

    Ok((pairs, row_genomes, col_genomes))
}

pub fn calc_allele_frequency(
    rec: &mut GenotypeRecords,
    num_genomes: usize,
) -> Result<AHashMap<(u32, u8), f64>> {
    if !rec.is_sorted_by_postion_allele_genome() {
        rec.sort_by_position_allele_genome()
            .change_context(GtencodeError::Input)?;
    }
    let posallele2freq = rec
        .records()
        .chunk_by(|a, b| a.get_pos_allele() < b.get_pos_allele())
        .map(|chunks| {
            (
                chunks[0].get_pos_allele(),
                chunks.len() as f64 / num_genomes as f64,
            )
        })
        .collect();

    Ok(posallele2freq)
}

pub fn calc_allele_count(rec: &mut GenotypeRecords) -> Result<AHashMap<u32, u32>> {
    rec.sort_by_position_genome_allele()
        .change_context(GtencodeError::Input)?;
    let mut count_map = AHashMap::<u32, u32>::new();
    for set in rec.records().linear_group_by_key(|x| x.get_position()) {
        let pos = set[0].get_position();
        let ac = set.len() as u32;
        count_map.insert(pos, ac);
    }
    Ok(count_map)
}
/// Calculate the base sum of relationship (assuming all genotypes are reference allele)
///
/// This is helpful as majority all the genotype are 0s (x=2 references).
/// We can iterative update the run sum of relatedships when we see a non-ref allelel:
/// subtracting a value assume ref allele and add a value consider non-ref
pub fn calc_base_relationship(freq_map: &AHashMap<u32, f64>) -> f64 {
    freq_map
        .values()
        .map(|p| (2.0 - 2.0 * p) * (2.0 - 2.0 * p) / 2.0 / p / (1.0 - p))
        .sum()
}

pub fn read_groups_file(
    p: impl AsRef<Path>,
    inds: &Individuals,
    level: Level,
) -> Result<HashMap<String, Vec<u32>>> {
    // parse group files
    let mut ind_id = 0;
    let mut grp_name = "";
    let mut group_map = HashMap::<String, Vec<u32>>::new();
    for line in std::fs::read_to_string(p.as_ref())
        .change_context(GtencodeError::Input)
        .attach("fail to read group files")?
        .trim()
        .split("\n")
    {
        ensure!(
            line.split("\t").count() == 2,
            GtencodeError::Input
                .into_report()
                .attach("number of columns is not two in groups file")
        );
        for (ifield, field) in line.split("\t").enumerate() {
            match ifield {
                0 => {
                    ind_id = if let Level::IndividualLevel = level {
                        *inds
                            .m()
                            .get(field)
                            .ok_or(GtencodeError::Input)
                            .attach("invalid sample name in groups file")?
                    } else {
                        field
                            .parse()
                            .change_context(GtencodeError::Input)
                            .attach("error in parsing genome id")?
                    };
                }
                1 => {
                    grp_name = field;
                }
                _ => {}
            }
        }
        if let Some(val) = group_map.get_mut(grp_name) {
            val.push(ind_id as u32);
        } else {
            group_map.insert(grp_name.to_owned(), vec![ind_id as u32]);
        }
    }
    Ok(group_map)
}

pub fn read_and_concat_rare_genotypes(
    records_paths: &[PathBuf],
    min_ac: u32,
    max_ac: u32,
) -> Result<(GenotypeRecords, Individuals)> {
    if records_paths.is_empty() {
        bail!(GtencodeError::Input
            .into_report()
            .attach("at least one record file should be specified"));
    }
    let mut records = GenotypeRecords::from_parquet_file(&records_paths[0])
        .change_context(GtencodeError::Input)?;
    filter_rv_by_frequency(&mut records, min_ac, max_ac)?;
    let ind_file = records_paths[0].with_extension("ind");
    let inds = Individuals::from_parquet_file(&ind_file)
        .change_context(GtencodeError::Input)
        .attach("fail to read individual file")?;
    for rec_file in records_paths.iter().skip(1) {
        let mut records_additional =
            GenotypeRecords::from_parquet_file(rec_file).change_context(GtencodeError::Input)?;
        filter_rv_by_frequency(&mut records_additional, min_ac, max_ac)?;
        let ind_file = rec_file.with_extension("ind");
        let inds_additional = Individuals::from_parquet_file(&ind_file)
            .change_context(GtencodeError::Input)
            .attach("fail to read individual file")?;
        ensure!(
            inds.v() == inds_additional.v(),
            GtencodeError::Input
                .into_report()
                .attach("individuals are not consistent across record files")
        );
        records.merge(records_additional);
    }
    Ok((records, inds))
}

pub fn filter_rv_by_frequency(recs: &mut GenotypeRecords, min_ac: u32, max_ac: u32) -> Result<()> {
    if !recs.is_sorted_by_postion_allele_genome() {
        recs.sort_by_position_allele_genome()
            .change_context(GtencodeError::Library)
            .attach("fail to sort rare genotype by position/allele/genome")?;
    }
    recs.records_mut()
        .linear_group_by_key_mut(|rec| rec.get_pos_allele())
        .filter(|blk| blk.len() < min_ac as usize || blk.len() > max_ac as usize)
        .for_each(|blk| blk.iter_mut().for_each(|rec| rec.set_sentinel()));

    recs.records_mut().retain(|rec| !rec.is_sentinel());
    Ok(())
}
