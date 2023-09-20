use ahash::AHashMap;
use ishare::genotype::rare::GenotypeRecords;
use itertools::Itertools;
use slice_group_by::GroupBy;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

pub fn file_to_u32_vec(p: impl AsRef<Path>) -> Vec<u32> {
    let mut s = String::new();
    let mut reader = std::fs::File::open(p).map(std::io::BufReader::new).unwrap();
    s.clear();
    use std::io::Read;
    reader.read_to_string(&mut s).unwrap();
    let mut v = s
        .trim()
        .split("\n")
        .map(|s| s.parse::<u32>().unwrap())
        .collect::<Vec<u32>>();
    assert!(v.len() != 0);
    v.sort();
    v
}

pub fn prep_pairs(
    records: &GenotypeRecords,
    genomes: &Option<Vec<u32>>,
    lists: &Vec<PathBuf>,
) -> (
    // pairs
    Vec<(u32, u32)>,
    // row genomes
    Vec<u32>,
    // col genomes
    Vec<u32>,
) {
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
            lst1.extend(file_to_u32_vec(&lists[0]));
            lst2.extend_from_slice(lst1.as_slice());
        }
        2 => {
            // inter-subset sharing
            lst1.extend(file_to_u32_vec(&lists[0]));
            lst2.extend(file_to_u32_vec(&lists[1]));

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

    let min_gid = records.records().first().unwrap().get_genome();
    let max_gid = records.records().last().unwrap().get_genome();
    let mut row_genomes = Vec::<u32>::new();
    let mut col_genomes = Vec::<u32>::new();
    let pairs: Vec<(u32, u32)> = match (genomes.len() > 0, lst1.len() > 0) {
        (true, false) => {
            row_genomes.extend(genomes.iter());
            col_genomes.extend(genomes.iter());
            genomes
                .iter()
                .map(|x| *x)
                .cartesian_product(genomes.iter().map(|x| *x))
                .filter(|(a, b)| *a > *b)
                .collect()
        }
        (false, true) => {
            row_genomes.extend(&lst1);
            col_genomes.extend(&lst2);

            if lists.len() == 2 {
                // two non-overallping list. no need to filter
                lst1.iter()
                    .map(|x| *x)
                    .cartesian_product(lst2.iter().map(|x| *x))
                    .collect()
            } else {
                // two identical list need to filtering
                lst1.iter()
                    .map(|x| *x)
                    .cartesian_product(lst2.iter().map(|x| *x))
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

    (pairs, row_genomes, col_genomes)
}

pub fn calc_allele_frequency(
    rec: &mut GenotypeRecords,
    num_hap: usize,
    num_sites: usize,
) -> AHashMap<u32, f64> {
    rec.sort_by_position();
    let mut freq_map = AHashMap::<u32, f64>::with_capacity(num_sites);
    let mut target_pos = u32::MAX;
    let mut cnt = 0u32;
    for r in rec.records() {
        let front_poistion = r.get_position();
        if front_poistion != target_pos {
            if target_pos != u32::MAX {
                // allele frequency
                let af = (cnt as f64) / (num_hap as f64);
                freq_map.insert(target_pos, af);
            }
            target_pos = front_poistion;
            cnt = 1;
        } else {
            cnt += 1;
        }
        if target_pos != u32::MAX {
            freq_map.insert(target_pos, (cnt as f64) / (num_hap as f64));
        }
    }
    freq_map
}

pub fn calc_allele_count(rec: &mut GenotypeRecords) -> AHashMap<u32, u32> {
    rec.sort_by_position();
    let mut count_map = AHashMap::<u32, u32>::new();
    for set in rec.records().linear_group_by_key(|x| x.get_position()) {
        let pos = set[0].get_position();
        let ac = set.len() as u32;
        count_map.insert(pos, ac);
    }
    count_map
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
