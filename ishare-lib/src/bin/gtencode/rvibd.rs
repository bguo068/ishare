use super::{GtencodeError, Result};
use error_stack::*;
use ishare::genome::Genome;
use ishare::share::mat::NamedMatrix;

use crate::utils::calc_allele_count;

use super::Commands;
use ahash::{AHashMap, AHashSet};
use ishare::container::intervaltree::IntervalTree;
use ishare::genotype::rare::GenotypeRecords;
use ishare::gmap::GeneticMap;
use ishare::indiv::Individuals;
use ishare::utils::path::from_prefix;
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

pub fn main_rvibd(args: &Commands) -> Result<()> {
    if let Commands::RvIBD {
        eibd,
        rec,
        samples_lst,
        genome_info,
        out_prefix,
        which,
    } = args
    {
        let valid_which = [0, 1, 2, 3];
        if !valid_which.iter().any(|x| x == which) {
            eprintln!("please specify the correct which!");
            std::process::exit(-1);
        }
        env_logger::Builder::new()
            .filter_level(LevelFilter::Info)
            .init();

        info!("read genome info"); // (for ibd) to compare with that from rare variants
        let (ginfo, gmap) = if genome_info
            .extension()
            .map(|e| e == "toml")
            .unwrap_or(false)
        {
            let ginfo =
                GenomeInfo::from_toml_file(genome_info).change_context(GtencodeError::Input)?;
            let gmap = GeneticMap::from_genome_info(&ginfo).change_context(GtencodeError::Input)?;
            (ginfo, gmap)
        } else {
            Genome::load_from_bincode_file(genome_info)
                .change_context(GtencodeError::Input)
                .attach("cannot load genome from binary flile")?
                .into_parts()
        };

        info!("read rv records");
        let mut records =
            GenotypeRecords::from_parquet_file(rec).change_context(GtencodeError::Input)?;

        info!("read individuals");
        let ind_file = rec.with_extension("ind");
        let inds =
            Individuals::from_parquet_file(&ind_file).change_context(GtencodeError::Input)?;

        info!("read ibd/rv samples");
        check_samples_orders(samples_lst, &inds)?;

        info!("read ibd");
        // ibd interval trees
        let mut ibd = read_ibdseg_vec(eibd).change_context(GtencodeError::Input)?;

        match *which {
            0 => position_scan(records, ibd, &ginfo, out_prefix)?,
            1 => pairwise_compare(records, ibd, &gmap, inds.v().len() * 2, out_prefix)?,
            2 => cmp_rv_and_ibd_similarity(
                &mut ibd,
                &mut records,
                inds.v().len() as u32 * 2,
                &gmap,
                out_prefix,
            )?,
            3 => {
                cmp_rv_and_ibd_len(&mut ibd, &mut records, &gmap, out_prefix)?;
            }
            _ => panic!("Not implemented"),
        }
    }
    Ok(())
}

fn check_samples_orders(samples_lst: &PathBuf, inds: &Individuals) -> Result<()> {
    let mut ibd_samples = vec![];
    std::fs::read_to_string(samples_lst)
        .change_context(GtencodeError::Input)?
        .trim()
        .split("\n")
        .for_each(|x| {
            ibd_samples.push(x.to_owned());
        });
    ensure!(
        inds.v() == &ibd_samples,
        GtencodeError::Input
            .into_report()
            .attach("sample order not matched")
    );
    Ok(())
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
) -> Result<()> {
    info!("remove multiallelic sites");
    records
        .filter_multi_allelic_site()
        .change_context(GtencodeError::Library)?;

    info!("sort rv records by position");
    records
        .sort_by_position()
        .change_context(GtencodeError::Library)?;
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
        let out =
            from_prefix(out_prefix.as_ref(), "rvibd.pos").change_context(GtencodeError::Output)?;
        // let out = format!("{}_rvibd.pos", out_prefix.as_ref().to_str().unwrap());
        File::create(&out)
            .map(BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .change_context(GtencodeError::Output)?
    };

    chunks.par_iter().try_for_each(|&(start, end)| {
        println!("=> {start} - {end}");
        position_scan_chunk(start, end, &records, ginfo, &tree, rwlock_file.clone())
    })?;
    Ok(())
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
) -> Result<()> {
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
    let mut file = rwlock_file
        .write()
        .map_err(|_| GtencodeError::Output)
        .attach("error with rwlock")?;
    for r in res {
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.chrid, r.chrpos, r.ac, r.within, r.between, r.out,
        )
        .change_context(GtencodeError::Output)?;
    }
    println!("<============ {start} - {end}");
    Ok(())
}

fn pairwise_compare(
    mut records: GenotypeRecords,
    mut ibd: Vec<IbdSeg>,
    gmap: &GeneticMap,
    nhap: usize,
    out_prefix: impl AsRef<Path>,
) -> Result<()> {
    info!("get allele count map");
    let ac_map = calc_allele_count(&mut records).change_context(GtencodeError::Library)?;

    info!("remove multiallelic sites");
    records
        .filter_multi_allelic_site()
        .change_context(GtencodeError::Library)?;

    info!("sort rv records by genome");
    records
        .sort_by_genome()
        .change_context(GtencodeError::Library)?;

    info!("sort ibd by genome pair");
    ibd.par_sort();

    let total_pair = nhap * (nhap - 1) / 2;

    let rwlock_file = {
        let out = out_prefix.as_ref().with_extension("rvibd.matched");
        File::create(&out)
            .map(BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .change_context(GtencodeError::Library)
    }?;

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
                gmap,
                &ac_map,
            )
        })
        .collect();

    let mut acc_counters = [0usize; 6];
    for c in res {
        let c = c?;
        acc_counters
            .iter_mut()
            .zip(c.iter())
            .for_each(|(acc, c)| *acc += *c);
    }

    let out = out_prefix.as_ref().with_extension("rvibd.pair");
    let mut file = File::create(&out)
        .map(BufWriter::new)
        .change_context(GtencodeError::Library)?;
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{}",
        acc_counters[0],
        acc_counters[1],
        acc_counters[2],
        acc_counters[3],
        acc_counters[4],
        acc_counters[5],
    )
    .change_context(GtencodeError::Library)?;
    println!("{acc_counters:?}");
    Ok(())
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
) -> Result<[usize; 6]> {
    info!("entering: {tri_idx1} - {tri_idx2}");
    let mut tree = IntervalTree::new(100);

    let mut counter = [0usize; 6];

    let pair_iter = (tri_idx1..tri_idx2).map(idx_to_i_j);

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
                    match (a, b, ele) {
                        (true, false, None) => counter[0] += 1,
                        (false, true, None) => counter[1] += 1,
                        (true, true, None) => counter[2] += 1,
                        (true, false, Some(_)) => counter[3] += 1,
                        (false, true, Some(_)) => counter[4] += 1,
                        (true, true, Some(ele)) => {
                            counter[5] += 1;
                            let rng = &ele.range;
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

    let mut file = rwlock_file
        .write()
        .map_err(|_| GtencodeError::Output)
        .attach("error accessing rwlock file")?;
    for r in v {
        file.write_all(&r.g1.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.g2.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.s.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.e.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.cm.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.p.to_le_bytes())
            .change_context(GtencodeError::Output)?;
        file.write_all(&r.ac.to_le_bytes())
            .change_context(GtencodeError::Output)?;
    }

    info!("leaving: {tri_idx1} - {tri_idx2}");
    Ok(counter)
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

/// The goal is generate a table, where rows are different IBD length bins,
/// and columns are for different AC values/value bins. The cells are the the
/// number of counts of allels shared that has AC and IBD length belongs each
/// categories.
fn cmp_rv_ac_and_ibd_len_chunks(
    ibd: &[IbdSeg],
    rvgt: &GenotypeRecords,
    rvgt2: &GenotypeRecords,
    pair_chunk: &[(u32, u32)],
    gmap: &GeneticMap,
    ibdlen_bins: &[f32],
    ac_bins: &[u32],
) -> Result<(Vec<u32>, NamedMatrix<u32>)> {
    // [0, 1, 2, 3], then we need to bins [-inf, 0), [0, 1), [1, 2), [2, 3) and [3, Inf)
    let n_ibdbins = ibdlen_bins.len() + 1;
    let n_ac = ac_bins.len() + 1;
    let mut mat = NamedMatrix::<u32>::new_from_shape(n_ibdbins as u32, n_ac as u32);
    let mut nonibd_ac = vec![0u32; n_ac];

    let ibdblk_iter = ibd.linear_group_by_key(|seg| seg.genome_id_pair());
    let merged_iter = pair_chunk
        .iter()
        .merge_join_by(ibdblk_iter, |&a, b| a.cmp(&b[0].genome_id_pair()));

    let mut tree = IntervalTree::new(100);

    for item in merged_iter {
        match item {
            itertools::EitherOrBoth::Both(pair, blk) => {
                let (genome1, genome2) = pair;
                // all IBD for a given pair is in the tree
                tree.clear_and_fill_with_iter(blk.iter().map(|seg| {
                    let cm = seg.get_seg_len_cm(gmap);
                    let ibdbin_idx = ibdlen_bins.partition_point(|x| *x <= cm);
                    (seg.s..seg.e, ibdbin_idx)
                }));

                for (pos, allele) in rvgt.iterate_rv_shared_for_genome_pair(*genome1, *genome2) {
                    let ac = rvgt2.get_allele_count_by_position_allele(pos, allele);
                    let ac_bin_ix = ac_bins.partition_point(|x| *x <= ac as u32);
                    if let Some(e) = tree.query_point(pos).next() {
                        let ibdlen_idx = e.value;
                        *mat.ref_mut_by_positions(ibdlen_idx as u32, ac_bin_ix as u32) += 1;
                    } else {
                        nonibd_ac[ac_bin_ix] += 1;
                    }
                }
            }
            itertools::EitherOrBoth::Left(pair) => {
                let (genome1, genome2) = pair;
                for (pos, allele) in rvgt.iterate_rv_shared_for_genome_pair(*genome1, *genome2) {
                    let ac = rvgt2.get_allele_count_by_position_allele(pos, allele);
                    let ac_bin_ix = ac_bins.partition_point(|x| *x <= ac as u32);
                    nonibd_ac[ac_bin_ix] += 1;
                }
            }
            _ => {}
        }
    }

    Ok((nonibd_ac, mat))
}

fn cmp_rv_and_ibd_similarity_chunks(
    ibd: &[IbdSeg],
    rvgt: &GenotypeRecords,
    pair_chunk: &[(u32, u32)],
    gmap: &GeneticMap,
    genome_span: f32,
    file: Arc<RwLock<BufWriter<File>>>,
) -> Result<()> {
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
    let mut file = file
        .write()
        .map_err(|_| GtencodeError::Output)
        .attach("error with rwlock")?;
    for r in v {
        writeln!(file, "{}\t{}\t{}", r.ibd_prop, r.cosine, r.jaccard)
            .change_context(GtencodeError::Output)?;
    }
    Ok(())
}
fn cmp_rv_and_ibd_similarity(
    ibd: &mut Vec<IbdSeg>,
    rvgt: &mut GenotypeRecords,
    nhap: u32,
    gmap: &GeneticMap,
    out_prefix: impl AsRef<Path>,
) -> Result<()> {
    let factor = 64;

    dbg!(nhap);

    // subsampling ibd and rv
    subsample(ibd, rvgt, factor);

    // sortting
    ibd.sort();
    rvgt.sort_by_genome()
        .change_context(GtencodeError::Library)?;

    // prepare output file
    let file = {
        let filename = out_prefix.as_ref().with_extension("pairwise_metrics.txt");
        // let filename = format!(
        //     "{}_pairwise_metrics.txt",
        //     out_prefix.as_ref().to_str().unwrap()
        // );
        File::create(&filename)
            .map(std::io::BufWriter::new)
            .map(RwLock::new)
            .map(Arc::new)
            .change_context(GtencodeError::Output)?
    };

    // some constant
    let genome_span = {
        let minmax = rvgt.records().iter().map(|x| x.get_position()).minmax();
        let (min, max) = minmax
            .into_option()
            .ok_or(GtencodeError::Library)
            .attach("can not find min/max position")?;
        gmap.get_cm_len(min, max)
    };

    // iterators
    // let npairs = nhap as usize * (nhap - 1) as usize;
    let pairs = (factor..nhap)
        .step_by(factor as usize)
        .flat_map(|i| (0..i).step_by(factor as usize).map(move |j| (i, j)))
        .filter(|(i, j)| (i % factor == 0) && (j % factor == 0) && (i > j))
        .collect_vec();

    pairs
        .chunks(200000)
        .par_bridge()
        .into_par_iter()
        .try_for_each(|pair_chunks| -> Result<()> {
            cmp_rv_and_ibd_similarity_chunks(
                ibd,
                rvgt,
                pair_chunks,
                gmap,
                genome_span,
                file.clone(),
            )?;
            Ok(())
        })?;

    eprintln!("genome_span = {genome_span}");
    Ok(())
}

fn cmp_rv_and_ibd_len(
    ibd: &mut [IbdSeg],
    rvgt: &mut GenotypeRecords,
    gmap: &GeneticMap,
    out_prefix: impl AsRef<Path>,
) -> Result<()> {
    // assumes: (1) IBD are sorted by genome pairs
    ensure!(
        ibd.is_sorted_by_key(|seg| (seg.genome_id_pair(), seg.coords())),
        GtencodeError::Library
            .into_report()
            .attach("IBD records are not sorted by genome pair and then by position")
    );

    // assumption 2: rvgt are sorted by genome pairs then by position and allele (use to find shared rv)
    ensure!(
        rvgt.is_sorted_by_genome()
            .change_context(GtencodeError::Library)
            .attach("error in checking if rvgt is properly sorted")?,
        GtencodeError::Library
            .into_report()
            .attach("rare genotype `rvgt` is not sorted by genome")
    );
    // assumption 3: rvgt2 are sorted by position and alleles (used to cout AC)
    let mut rvgt2 = rvgt.clone();
    rvgt2
        .sort_by_position()
        .change_context(GtencodeError::Library)
        .attach("fail to sort rare genotype by position")?;

    // get total number of haplotype/genomes
    let nhap = rvgt
        .records()
        .last()
        .map(|rec| rec.get_genome())
        .unwrap_or(0)
        + 1;

    // sortting
    rvgt.sort_by_genome()
        .change_context(GtencodeError::Library)?;

    let pairs = (1..nhap)
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .filter(|(i, j)| i > j)
        .collect_vec();
    let ibdlens = [0.0f32, 2.0, 3.0, 4.0, 6.0, 10.0, 18.0, 30.0];
    let ac_bins = (0..51).collect_vec();

    let mut res_vec: Vec<_> = pairs
        .chunks(2000)
        .par_bridge()
        .into_par_iter()
        .flat_map(|pair_chunks| -> Result<(Vec<u32>, NamedMatrix<u32>)> {
            cmp_rv_ac_and_ibd_len_chunks(ibd, rvgt, &rvgt2, pair_chunks, gmap, &ibdlens, &ac_bins)
        })
        .collect::<Vec<_>>();

    ensure!(
        pairs.len().div_ceil(2000) == res_vec.len(),
        GtencodeError::Library
            .into_report()
            .attach("not all chunks successfully completed rvs-ac-vs-IBD-len comparison analysis")
    );

    // combine data to the first element
    let (first, rest) = res_vec
        .split_first_mut()
        .ok_or(GtencodeError::Library)
        .attach("Empty result vector")?;
    let (nonibd_vec, ibd_mat) = first;
    for (nonibd_vec2, ibd_mat2) in rest {
        nonibd_vec
            .iter_mut()
            .zip(nonibd_vec2.iter())
            .for_each(|(a, b)| *a += *b);
        ibd_mat
            .get_data_slice_mut()
            .iter_mut()
            .zip(ibd_mat2.get_data_slice().iter())
            .for_each(|(a, b)| *a += *b);
    }
    // prepare output file
    let mut file = {
        let filename = out_prefix.as_ref().with_extension("pairwise_metrics.txt");
        File::create(&filename)
            .map(std::io::BufWriter::new)
            .change_context(GtencodeError::Output)?
    };

    // write header
    write!(file, "ibdlen_bin").change_context(GtencodeError::Output)?;
    for bin in &ac_bins {
        write!(file, "\tAC<{bin}").change_context(GtencodeError::Output)?;
    }
    writeln!(file, "\tAC<Inf").change_context(GtencodeError::Output)?;

    // write the non-ibd counts (first row)
    // write index first, then counts for ac bins
    write!(file, "NonIBD").change_context(GtencodeError::Output)?;
    for cnt in nonibd_vec.iter() {
        write!(file, "\t{cnt}").change_context(GtencodeError::Output)?;
    }
    writeln!(file).change_context(GtencodeError::Output)?; // add a newline

    // write IBD counts (non-first row)
    for (lenbin, ibd_vec) in ibdlens
        .iter()
        .zip(ibd_mat.get_data_slice().chunks(ac_bins.len() + 1))
    {
        write!(file, "IBD<{lenbin}").change_context(GtencodeError::Output)?;
        for cnt in ibd_vec.iter() {
            write!(file, "\t{cnt}").change_context(GtencodeError::Output)?;
        }
        writeln!(file).change_context(GtencodeError::Output)?; // add a newline
    }

    Ok(())
}
