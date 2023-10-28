use crate::gtencode_imp::utils::calc_allele_count;

use super::super::Commands;
use ahash::{AHashMap, AHashSet};
use ishare::container::intervaltree::IntervalTree;
use ishare::genotype::rare::GenotypeRecords;
use ishare::gmap::GeneticMap;
use ishare::indiv::Individuals;
use ishare::{genome::GenomeInfo, share::ibd::ibdseg::*};
use itertools::Itertools;
use log::{info, LevelFilter};
use rayon::prelude::*;
use slice_group_by::*;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};

trait GenomeIDPair {
    fn genome_id_pair(&self) -> (u32, u32);
}
impl GenomeIDPair for IbdSeg {
    fn genome_id_pair(&self) -> (u32, u32) {
        let (sid1, hid1, sid2, hid2) = self.haplotype_pair();
        let gid1 = (sid1 << 1) + hid1 as u32;
        let gid2 = (sid2 << 1) + hid2 as u32;
        (gid1, gid2)
    }
}

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
        let valid_which = [0, 1, 2];
        if !valid_which.iter().any(|x| x == which) {
            eprintln!("please specify the correct which!");
            std::process::exit(-1);
        }
        env_logger::Builder::new()
            .filter_level(LevelFilter::Info)
            .init();

        info!("read genome info"); // (for ibd) to compare with that from rare variants
        let ginfo = GenomeInfo::from_toml_file(genome_info);
        let gmap = GeneticMap::from_genome_info(&ginfo);

        info!("read rv records");
        let mut records = GenotypeRecords::from_parquet_file(rec);

        info!("read individuals");
        let ind_file = rec.with_extension("ind");
        let inds = Individuals::from_parquet_file(&ind_file);

        info!("read ibd/rv samples");
        check_samples_orders(samples_lst, &inds);

        info!("read ibd");
        // ibd interval trees
        let mut ibd = read_ibdseg_vec(eibd);

        match *which {
            0 => position_scan(records, ibd, &ginfo, out_prefix),
            1 => pairwise_compare(records, ibd, &gmap, inds.v().len() * 2, out_prefix),
            2 => cmp_rv_and_ibd_similarity(
                &mut ibd,
                &mut records,
                inds.v().len() as u32 * 2,
                &gmap,
                out_prefix,
            ),
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

    let it = ibd.into_iter().map(|seg| {
        let rng = seg.s..seg.e;
        let genome_pair = seg.genome_id_pair();
        (rng, genome_pair)
    });
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

    // prepare rwlock file

    let rwlock_file = {
        let out = format!("{}_rvibd.pos", out_prefix.as_ref().to_str().unwrap());
        File::create(&out)
            .map(BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .unwrap()
    };

    chunks.par_iter().for_each(|&(start, end)| {
        println!("=> {start} - {end}");
        position_scan_chunk(start, end, &records, &ginfo, &tree, rwlock_file.clone())
    });
}

/// Note on the tree argument: the range are IBD segment start/end (genome-wide
/// coordinates) the value is the genome-pair ids
fn position_scan_chunk(
    start: usize,
    end: usize,
    records: &GenotypeRecords,
    ginfo: &GenomeInfo,
    tree: &IntervalTree<u32, (u32, u32)>,
    rwlock_file: Arc<RwLock<BufWriter<File>>>,
) {
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
    let mut file = rwlock_file.write().unwrap();
    for r in res {
        write!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.chrid, r.chrpos, r.ac, r.within, r.between, r.out,
        )
        .unwrap();
    }
    println!("<============ {start} - {end}");
}

fn pairwise_compare(
    mut records: GenotypeRecords,
    mut ibd: Vec<IbdSeg>,
    gmap: &GeneticMap,
    nhap: usize,
    out_prefix: impl AsRef<Path>,
) {
    info!("get allele count map");
    let ac_map = calc_allele_count(&mut records);

    info!("remove multiallelic sites");
    records.filter_multi_allelic_site();

    info!("sort rv records by genome");
    records.sort_by_genome();

    info!("sort ibd by genome pair");
    ibd.par_sort();

    let total_pair = nhap * (nhap - 1) / 2;

    let rwlock_file = {
        let out = format!("{}_rvibd.matched", out_prefix.as_ref().to_str().unwrap());
        File::create(&out)
            .map(BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .unwrap()
    };

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
        .map(|(tri_idx1, tri_idx2)| {
            pairwise_compare_chunk(
                &records,
                &ibd,
                *tri_idx1,
                *tri_idx2,
                rwlock_file.clone(),
                &gmap,
                &ac_map,
            )
        })
        .collect();

    let mut acc_counters = [0usize; 6];
    for c in res {
        acc_counters
            .iter_mut()
            .zip(c.iter())
            .for_each(|(acc, c)| *acc += *c);
    }

    let out = format!("{}_rvibd.pair", out_prefix.as_ref().to_str().unwrap());
    let mut file = File::create(&out).map(BufWriter::new).unwrap();
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

struct MatchedRvIBD {
    g1: u32,
    g2: u32,
    s: u32,
    e: u32,
    cm: f32,
    p: u32,
    ac: u32,
}

fn pairwise_compare_chunk(
    records: &GenotypeRecords,
    ibd: &Vec<IbdSeg>,
    tri_idx1: usize,
    tri_idx2: usize,
    rwlock_file: Arc<RwLock<BufWriter<File>>>,
    gmap: &GeneticMap,
    ac_map: &AHashMap<u32, u32>,
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
            .partition_point(|x| x.genome_id_pair() < (i1, j1));
        let b = ibd
            .as_slice()
            .partition_point(|x| x.genome_id_pair() < (i2, j2));
        &ibd[a..b]
    };

    let ibdblk_iter = ibdslice.linear_group_by_key(|seg| seg.genome_id_pair());

    let merged_iter = pair_iter.merge_join_by(ibdblk_iter, |a, b| a.cmp(&b[0].genome_id_pair()));

    let mut v = Vec::<MatchedRvIBD>::new();
    for item in merged_iter {
        match item {
            itertools::EitherOrBoth::Both(pair, blk) => {
                // add ibd segs for a given genome pair to an interval tree
                let it = blk.iter().map(|seg| (seg.s..seg.e, ()));
                tree.clear_and_fill_with_iter(it);

                let (i, j) = pair;
                for (pos, a, b) in records.iter_genome_pair_genotypes(i, j) {
                    // whehter genome a has rv
                    let a = a.is_some();
                    // whether genome b has rv
                    let b = b.is_some();
                    // whether genome a and b share ibd over this rv
                    let ele = tree.query_point(pos).next();
                    let c = ele.is_some();
                    match (a, b, c) {
                        (true, false, false) => counter[0] += 1,
                        (false, true, false) => counter[1] += 1,
                        (true, true, false) => counter[2] += 1,
                        (true, false, true) => counter[3] += 1,
                        (false, true, true) => counter[4] += 1,
                        (true, true, true) => {
                            counter[5] += 1;
                            let rng = &ele.unwrap().range;
                            let matched = MatchedRvIBD {
                                g1: i,
                                g2: j,
                                s: rng.start,
                                e: rng.end,
                                p: pos,
                                cm: gmap.get_cm_len(rng.start, rng.end),
                                ac: ac_map[&pos],
                            };
                            v.push(matched);
                        }
                        _ => {}
                    }
                }
            }
            itertools::EitherOrBoth::Left(pair) => {
                // println!("left");
                let (i, j) = pair;
                for (_pos, a, b) in records.iter_genome_pair_genotypes(i, j) {
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
                // panic!("not possible");
            }
        }
    }

    let mut file = rwlock_file.write().unwrap();
    for r in v {
        file.write_all(&r.g1.to_le_bytes()).unwrap();
        file.write_all(&r.g2.to_le_bytes()).unwrap();
        file.write_all(&r.s.to_le_bytes()).unwrap();
        file.write_all(&r.e.to_le_bytes()).unwrap();
        file.write_all(&r.cm.to_le_bytes()).unwrap();
        file.write_all(&r.p.to_le_bytes()).unwrap();
        file.write_all(&r.ac.to_le_bytes()).unwrap();
    }

    info!("leaving: {} - {}", tri_idx1, tri_idx2);
    counter
}

fn subsample(ibd: &mut Vec<IbdSeg>, rvgt: &mut GenotypeRecords, factor: u32) {
    // mark for removal only keep those with genome ids are multiples of factor
    ibd.iter_mut()
        .filter(|seg| {
            let (i, j) = seg.genome_id_pair();
            (i % factor > 0) || (j % factor > 0)
        })
        .for_each(|seg| {
            seg.i = u32::MAX;
        });
    ibd.retain(|seg| seg.i != u32::MAX);

    // dbg!(&ibd[0..50].iter().map(|x| x.genome_id_pair()).collect_vec());
    rvgt.records_mut()
        .iter_mut()
        .filter(|r| r.get_genome() % factor != 0)
        .for_each(|r| r.set_sentinel());
    rvgt.records_mut().retain(|r| !r.is_sentinel());
    // dbg!(&rvgt.records().iter().map(|x| x.get_genome()).collect_vec());
}

struct MetricRecord {
    ibd_prop: f32,
    cosine: f32,
    jaccard: f32,
}

fn cmp_rv_and_ibd_similarity_chunks(
    ibd: &Vec<IbdSeg>,
    rvgt: &GenotypeRecords,
    pair_chunk: &[(u32, u32)],
    gmap: &GeneticMap,
    genome_span: f32,
    file: Arc<RwLock<BufWriter<File>>>,
) {
    let mut tree = IntervalTree::new(100);
    let ibdblk_iter = ibd.linear_group_by_key(|seg| seg.genome_id_pair());
    let merged_iter = pair_chunk
        .iter()
        .merge_join_by(ibdblk_iter, |&a, b| a.cmp(&b[0].genome_id_pair()));

    let mut v = vec![];
    // go over each pair of genomes
    for item in merged_iter {
        // dbg!(&item);
        let mut totibd = 0.0;
        let mut a = 0; // number of rare variants that g1 has
        let mut b = 0; // number of rare variatnts that g2 has
        let mut ab = 0; // number of rare variants that both g1 and g2 has

        let mut process_rare_variant_for_pair = |pair| {
            let (i, j) = pair;
            for (_pos, left, right) in rvgt.iter_genome_pair_genotypes(i, j) {
                match (left, right) {
                    (Some(_), Some(_)) => {
                        ab += 1;
                        a += 1;
                        b += 1;
                    }
                    (None, Some(_)) => b += 1,
                    (Some(_), None) => a += 1,
                    _ => {}
                }
            }
        };

        match item {
            itertools::EitherOrBoth::Both(pair, blk) => {
                // add ibd seg to for a given genome pair to an interval tree
                let it = blk.iter().map(|seg| (seg.s..seg.e, ()));
                tree.clear_and_fill_with_iter(it);

                println!("both");
                // calculate total ibd
                blk.iter()
                    .for_each(|seg| totibd += seg.get_seg_len_cm(gmap));

                process_rare_variant_for_pair(*pair);
            }
            itertools::EitherOrBoth::Left(pair) => {
                process_rare_variant_for_pair(*pair);
            }
            _ => {}
        }

        let ibd_prop = totibd / genome_span;
        let cosine = ab as f32 / ((a * b) as f32).sqrt();
        let jaccard = ab as f32 / (a + b - ab) as f32;

        v.push(MetricRecord {
            ibd_prop,
            cosine,
            jaccard,
        });
    }
    // write to file
    let mut file = file.write().unwrap();
    for r in v {
        write!(file, "{}\t{}\t{}\n", r.ibd_prop, r.cosine, r.jaccard).unwrap();
    }
}
fn cmp_rv_and_ibd_similarity(
    ibd: &mut Vec<IbdSeg>,
    rvgt: &mut GenotypeRecords,
    nhap: u32,
    gmap: &GeneticMap,
    out_prefix: impl AsRef<Path>,
) {
    let factor = 64;

    dbg!(nhap);

    // subsampling ibd and rv
    subsample(ibd, rvgt, factor);

    // sortting
    ibd.sort();
    rvgt.sort_by_genome();

    // prepare output file
    let file = {
        let filename = format!(
            "{}_pairwise_metrics.txt",
            out_prefix.as_ref().to_str().unwrap()
        );
        File::create(&filename)
            .map(std::io::BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .unwrap()
    };

    // some constant
    let genome_span = {
        let minmax = rvgt.records().iter().map(|x| x.get_position()).minmax();
        let (min, max) = minmax.into_option().unwrap();
        gmap.get_cm_len(min, max)
    };

    // iterators
    // let npairs = nhap as usize * (nhap - 1) as usize;
    let pairs = (factor..nhap)
        .into_iter()
        .step_by(factor as usize)
        .flat_map(|i| {
            (0..i)
                .into_iter()
                .step_by(factor as usize)
                .map(move |j| (i, j))
        })
        .filter(|(i, j)| (i % factor == 0) && (j % factor == 0) && (i > j))
        .collect_vec();

    pairs
        .chunks(200000)
        .par_bridge()
        .into_par_iter()
        .for_each(|pair_chunks| {
            cmp_rv_and_ibd_similarity_chunks(
                &ibd,
                &rvgt,
                pair_chunks,
                gmap,
                genome_span,
                file.clone(),
            );
        });

    eprintln!("genome_span = {genome_span}");
}
