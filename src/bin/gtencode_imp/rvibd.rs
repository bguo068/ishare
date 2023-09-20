use super::super::Commands;
use ahash::AHashSet;
use ishare::container::intervaltree::IntervalTree;
use ishare::genotype::rare::GenotypeRecords;
use ishare::indiv::Individuals;
use ishare::{genome::GenomeInfo, share::ibd::ibdseg::*};
use itertools::Itertools;
use log::{info, LevelFilter};
use rayon::prelude::*;
use slice_group_by::*;
use std::path::{Path, PathBuf};

pub fn main_rvibd(args: &Commands) {
    if let Commands::RvIBD {
        eibd,
        rec,
        samples_lst,
        genome_info,
        out_prefix,
        which,
    } = args
    {
        let valid_which = [0, 1];
        if !valid_which.iter().any(|x| x == which) {
            eprintln!("please specify the correct which!");
            std::process::exit(-1);
        }
        env_logger::Builder::new()
            .filter_level(LevelFilter::Info)
            .init();

        info!("read genome info"); // (for ibd) to compare with that from rare variants
        let ginfo = GenomeInfo::from_toml_file(genome_info);

        info!("read rv records");
        let records = GenotypeRecords::from_parquet_file(rec);

        info!("read individuals");
        let ind_file = rec.with_extension("ind");
        let inds = Individuals::from_parquet_file(&ind_file);

        info!("read ibd/rv samples");
        check_samples_orders(samples_lst, &inds);

        info!("read ibd");
        // ibd interval trees
        let ibd = read_ibdseg_vec(eibd);

        match *which {
            0 => position_scan(records, ibd, &ginfo, out_prefix),
            1 => pairwise_compare(records, ibd, &ginfo, inds.v().len() * 2, out_prefix),
            _ => panic!("Not implemented"),
        }
    }
}

fn check_samples_orders(samples_lst: &PathBuf, inds: &Individuals) {
    let mut ibd_samples = vec![];
    std::fs::read_to_string(samples_lst)
        .expect("can not open sample list file")
        .trim()
        .split("\n")
        .for_each(|x| {
            ibd_samples.push(x.to_owned());
        });
    assert_eq!(inds.v(), &ibd_samples);
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct PositionScanRecord {
    chrid: u32,
    chrpos: u32,
    ac: u32,
    within: usize,
    between: usize,
    out: usize,
}

fn position_scan(
    mut records: GenotypeRecords,
    ibd: Vec<IbdSeg>,
    ginfo: &GenomeInfo,
    out_prefix: impl AsRef<Path>,
) {
    info!("remove multiallelic sites");
    records.filter_multi_allelic_site();

    info!("sort rv records by position");
    records.sort_by_position();
    info!("build ibd interval tress");
    let it = ibd.into_iter().map(|seg| (seg.s..seg.e, (seg.i, seg.j)));
    let tree = IntervalTree::from_iter(it);

    dbg!(tree.iter().count());
    dbg!(records.records().len());

    info!("iter each position as a chunk");

    let chunks = {
        let step = records.records().len() / 1000;
        let mut start = 0;
        let mut end = 0;
        let mut v = vec![];

        records
            .records()
            .linear_group_by_key(|x| x.get_position())
            .for_each(|blk| {
                end += blk.len();
                if end - start > step {
                    v.push((start, end));
                    start = end;
                }
            });
        if start != end {
            v.push((start, end));
        }
        v
    };

    let res: Vec<_> = chunks
        .par_iter()
        .map(|&(start, end)| {
            println!("=> {start} - {end}");
            position_scan_chunk(start, end, &records, &ginfo, &tree)
        })
        .collect();
    let mut v: Vec<_> = res.into_iter().flat_map(|v| v.into_iter()).collect();
    v.sort();

    let out = format!("{}_rvibd.psc", out_prefix.as_ref().to_str().unwrap());
    let mut file = std::fs::File::create(&out)
        .map(std::io::BufWriter::new)
        .unwrap();
    use std::io::Write;
    for r in v {
        write!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.chrid, r.chrpos, r.ac, r.within, r.between, r.out,
        )
        .unwrap();
    }
}

fn position_scan_chunk(
    start: usize,
    end: usize,
    records: &GenotypeRecords,
    ginfo: &GenomeInfo,
    tree: &IntervalTree<u32, (u32, u32)>,
) -> Vec<PositionScanRecord> {
    let mut genomes_with_rv = AHashSet::<u32>::new();

    let mut res = Vec::new();
    for rv in records.records()[start..end]
        .linear_group_by_key(|x| x.get_position())
        .filter(|blk| blk.len() > 1)
    {
        genomes_with_rv.clear();
        genomes_with_rv.extend(rv.iter().map(|x| x.get_genome()));
        let pos = rv[0].get_position();
        let (chrid, _chrname, chr_pos) = ginfo.to_chr_pos(pos);
        let ac = genomes_with_rv.len();

        let mut within = 0;
        let mut between = 0;
        let mut out = 0;

        for elm in tree.query_point(pos) {
            let (i, j) = elm.value;
            match (genomes_with_rv.contains(&i), genomes_with_rv.contains(&j)) {
                (true, true) => within += 1,
                (false, false) => out += 1,
                _ => between += 1,
            }
        }

        res.push(PositionScanRecord {
            chrid: chrid as u32,
            chrpos: chr_pos,
            ac: ac as u32,
            within,
            between,
            out,
        });
    }
    println!("<============ {start} - {end}");

    res
}

fn pairwise_compare(
    mut records: GenotypeRecords,
    mut ibd: Vec<IbdSeg>,
    _ginfo: &GenomeInfo,
    nhap: usize,
    out_prefix: impl AsRef<Path>,
) {
    info!("sort rv records by genome");
    records.sort_by_genome();

    info!("sort ibd by genome pair");
    ibd.par_sort();

    let total_pair = nhap * (nhap - 1) / 2;

    // split into ~ 1000 chunks
    let chunks = {
        let mut chunks = vec![];
        let step = total_pair / 1000;
        let mut s = 0;
        let mut e;
        while s < total_pair {
            e = s + step;
            if e > total_pair {
                e = total_pair;
            }
            chunks.push((s, e));
            s = e;
        }
        chunks
    };

    let res: Vec<_> = chunks
        .par_iter()
        .map(|(tri_idx1, tri_idx2)| pairwise_compare_chunk(&records, &ibd, *tri_idx1, *tri_idx2))
        .collect();

    let mut acc_counters = [0usize; 6];
    for c in res {
        acc_counters
            .iter_mut()
            .zip(c.iter())
            .for_each(|(acc, c)| *acc += *c);
    }

    let out = format!("{}_rvibd.pwc", out_prefix.as_ref().to_str().unwrap());
    let mut file = std::fs::File::create(&out)
        .map(std::io::BufWriter::new)
        .unwrap();
    use std::io::Write;
    write!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{}\n",
        acc_counters[0],
        acc_counters[1],
        acc_counters[2],
        acc_counters[3],
        acc_counters[4],
        acc_counters[5],
    )
    .unwrap();
    println!("{:?}", acc_counters);
}

fn idx_to_i_j(idx: usize) -> (u32, u32) {
    let row = ((2.0 * idx as f64 + 0.25).sqrt() + 0.5) as u32;
    // overflow if u32_t * u32_t. Need to convert to usize
    let tmp = row as usize;
    let col = (idx - tmp * (tmp - 1) / 2) as u32;
    (row, col)
}

#[test]
fn test_idx_to_i_j() {
    assert_eq!(idx_to_i_j(0), (1, 0));
    assert_eq!(idx_to_i_j(1), (2, 0));
    assert_eq!(idx_to_i_j(2), (2, 1));
    assert_eq!(idx_to_i_j(3), (3, 0));
    assert_eq!(idx_to_i_j(4), (3, 1));
    assert_eq!(idx_to_i_j(5), (3, 2));
}
fn pairwise_compare_chunk(
    records: &GenotypeRecords,
    ibd: &Vec<IbdSeg>,
    tri_idx1: usize,
    tri_idx2: usize,
) -> [usize; 6] {
    info!("entering: {} - {}", tri_idx1, tri_idx2);
    let mut tree = IntervalTree::new(100);

    let mut counter = [0usize; 6];

    let pair_iter = (tri_idx1..tri_idx2).into_iter().map(|idx| idx_to_i_j(idx));

    // find ibd slice
    let ibdslice = {
        let (i1, j1) = idx_to_i_j(tri_idx1);
        let (i2, j2) = idx_to_i_j(tri_idx2);
        let a = ibd
            .as_slice()
            .partition_point(|x| x.haplotype_pair_int() < (i1, j1));
        let b = ibd
            .as_slice()
            .partition_point(|x| x.haplotype_pair_int() < (i2, j2));
        &ibd[a..b]
    };

    let ibdblk_iter = ibdslice.linear_group_by_key(|seg| seg.haplotype_pair_int());

    let merged_iter = pair_iter.merge_join_by(ibdblk_iter, |a, b| a.cmp(&(b[0].s, b[0].e)));

    for item in merged_iter {
        match item {
            itertools::EitherOrBoth::Both(_pair, blk) => {
                // add ibd seg to tree
                let it = blk.iter().map(|seg| (seg.s..seg.e, ()));
                tree.clear_and_fill_with_iter(it);

                let (i, j) = blk[0].haplotype_pair_int();
                for (pos, a, b) in records.iter_genome_pair_genotypes(i, j) {
                    let a = a.is_some();
                    let b = b.is_some();
                    let c = tree.query_point(pos).next().is_some();
                    match (a, b, c) {
                        (true, false, false) => counter[0] += 1,
                        (false, true, false) => counter[1] += 1,
                        (true, true, false) => counter[2] += 1,
                        (true, false, true) => counter[3] += 1,
                        (false, true, true) => counter[4] += 1,
                        (true, true, true) => counter[5] += 1,
                        _ => {}
                    }
                }
            }
            itertools::EitherOrBoth::Left(pair) => {
                for (_pos, a, b) in records.iter_genome_pair_genotypes(pair.0, pair.1) {
                    let a = a.is_some();
                    let b = b.is_some();
                    let c = false;
                    match (a, b, c) {
                        (true, false, false) => counter[0] += 1,
                        (false, true, false) => counter[1] += 1,
                        (true, true, false) => counter[2] += 1,
                        _ => {}
                    }
                }
            }
            itertools::EitherOrBoth::Right(_blk) => {
                panic!("not possible");
            }
        }
    }
    info!("leaving: {} - {}", tri_idx1, tri_idx2);
    counter
}
