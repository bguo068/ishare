use std::num::ParseFloatError;

use super::utils::*;
use ishare::{
    genome,
    gmap::{self, GeneticMap},
    indiv::*,
    share::ibd::{
        ibdseg::IbdSeg,
        ibdset::*,
        overlap::{self, write_per_winddow_overlap_res, IbdOverlapResult},
    },
};

use itertools::Itertools;
use slice_group_by::GroupBy;
use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Indiv { source: ishare::indiv::Error },
    #[snafu(transparent)]
    Genome { source: ishare::genome::Error },
    #[snafu(transparent)]
    Gmap { source: ishare::gmap::Error },
    #[snafu(transparent)]
    Ibd { source: ishare::share::ibd::Error },
    #[snafu(transparent)]
    IbdutilsUtil { source: super::utils::Error },
    // local
    #[snafu(transparent)]
    ParseFloat { source: ParseFloatError },
    #[snafu(transparent)]
    StdIo { source: std::io::Error },
}
type Result<T> = std::result::Result<T, Error>;

use super::super::Commands;
pub fn main_compare(args: &Commands) -> Result<()> {
    if let Commands::Compare {
        genome_info,
        sample_lst1,
        sample_lst2,
        fmt1,
        fmt2,
        ibd1_dir,
        ibd2_dir,
        min_cm,
        length_bin_starts,
        window_size_bp,
        use_hap_overlap,
        use_hap_totibd,
        suppress_total_ibd_calculation,
        write_details,
        out,
    } = args
    {
        // files
        let ginfo = genome::GenomeInfo::from_toml_file(genome_info)?;
        let gmap = gmap::GeneticMap::from_genome_info(&ginfo)?;

        let (inds1, inds1_opt) = Individuals::from_txt_file(sample_lst1)?;
        let (inds2, inds2_opt) = Individuals::from_txt_file(sample_lst2)?;
        let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds1);
        let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds2);

        for (((ibd, fmt), dir), inds_opt) in [&mut ibd1, &mut ibd2]
            .iter_mut()
            .zip([fmt1, fmt2])
            .zip([ibd1_dir, ibd2_dir])
            .zip([&inds1_opt, &inds2_opt])
        {
            if fmt.as_str() == "hapibd" {
                ibd.read_hapibd_dir(dir)?;
                match *use_hap_overlap {
                    true => ibd.sort_by_haplotypes(),
                    false => ibd.sort_by_samples(),
                }
                ibd.infer_ploidy();
            } else if fmt.as_str() == "tskibd" {
                ibd.read_tskibd_dir(dir)?;
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else if fmt.as_str() == "hmmibd" {
                ibd.read_hmmibd_dir(dir)?;
                ibd.sort_by_haplotypes();
                ibd.infer_ploidy();
            } else {
                panic!("format {fmt} is not supported.");
            }

            // remove short IBD segment to prevent unfair genome-wide overlapping analysis
            // for instance, original hmmIBD output all short IBD segments without filtering;
            // while other ibd caller such hap-IBD by default only output cm >= 2.0.
            // Directly filtering IBD segments by cm length take additonal bp to cM mapping.
            // Adding it here seems to be convenient
            ibd.filter_segments_by_min_cm(*min_cm);

            match inds_opt.as_ref() {
                Some((converter, ind_, PloidConvertDirection::Diploid2Haploid)) => {
                    ibd.covert_to_haploid(ind_, converter);
                }
                Some((converter, ind_, PloidConvertDirection::Haploid2Diploid)) => {
                    ibd.covert_to_het_diploid(ind_, converter)?;
                }
                None => {}
            }
        }

        assert_eq!(ibd1.get_inds().v(), ibd2.get_inds().v());

        {
            // convert from string to vector of float
            let mut v: Vec<f32> = vec![];
            length_bin_starts.split(",").try_for_each(|s| {
                v.push(s.parse()?);
                Ok::<(), Error>(())
            })?;
            let mut length_bin_starts = v;
            // sort
            length_bin_starts.sort_by(|a, b| a.total_cmp(b));

            let ignore_hap = !(*use_hap_overlap);

            // 1. overlapping analysis
            let prefix_for_details = match *write_details {
                true => Some(out),
                false => None,
            };
            let oa = overlap::IbdOverlapAnalyzer::new(&ibd1, &ibd2, ignore_hap, prefix_for_details);

            // 1.1 overlapping analysis all together
            let res = oa.analzyze(Some(length_bin_starts.as_slice()), None)?;
            res.to_csv(out)?;

            // 1.2 overlapping analysis per window of each chromosome
            if let Some(window_size_bp) = window_size_bp {
                let nchrom = ginfo.gwstarts.len();
                let mut windows = Vec::<(u32, u32)>::new();
                let mut ov_res_windows = Vec::<IbdOverlapResult>::new();
                for idx in 0..nchrom {
                    let chr_gw_start = ginfo.gwstarts[idx];
                    let chr_gw_end = if idx == nchrom - 1 {
                        ginfo.chromsize[nchrom - 1] + ginfo.gwstarts[nchrom - 1]
                    } else {
                        ginfo.gwstarts[idx + 1]
                    };
                    let mut winstart = chr_gw_start;
                    while winstart < chr_gw_end {
                        let mut winend = winstart + window_size_bp - 1;
                        if winend > chr_gw_end {
                            winend = chr_gw_end - 1;
                        }
                        windows.push((winstart, winend));
                        let ov_res = oa.analzyze(
                            Some(length_bin_starts.as_slice()),
                            Some((winstart, winend)),
                        )?;
                        ov_res_windows.push(ov_res);
                        winstart += window_size_bp;
                    }
                }
                // write to file
                let out_win_ov = out.with_extension("winovcsv");
                write_per_winddow_overlap_res(
                    &windows,
                    &ov_res_windows,
                    &ginfo,
                    &gmap,
                    &out_win_ov,
                )?;
            }
        }
        if !suppress_total_ibd_calculation {
            // 2. total IBD analysis
            let ignore_hap = !(*use_hap_totibd);
            let mut total1_vec = vec![];
            let mut total2_vec = vec![];

            // re-sort ibd
            match ignore_hap {
                true => {
                    for ibd in [&mut ibd1, &mut ibd2] {
                        if !ibd.is_sorted_by_samples() {
                            ibd.sort_by_samples();
                        }
                        ibd.infer_ploidy();
                    }
                    let it1 = ibd1
                        .as_slice()
                        .linear_group_by_key(|seg| seg.individual_pair());
                    let it2 = ibd2
                        .as_slice()
                        .linear_group_by_key(|seg| seg.individual_pair());
                    for e in it1.merge_join_by(it2, |blk1, blk2| {
                        blk1[0].individual_pair().cmp(&blk2[0].individual_pair())
                    }) {
                        // individual pairs. segments may overlap
                        match e {
                            itertools::EitherOrBoth::Left(blk1) => {
                                let total1: f32 = get_sample_pair_total_ibd(blk1, &gmap);
                                total1_vec.push(total1);
                                total2_vec.push(0.0f32);
                            }
                            itertools::EitherOrBoth::Right(blk2) => {
                                let total2: f32 = get_sample_pair_total_ibd(blk2, &gmap);
                                total2_vec.push(total2);
                                total1_vec.push(0.0f32);
                            }
                            itertools::EitherOrBoth::Both(blk1, blk2) => {
                                let total1: f32 = get_sample_pair_total_ibd(blk1, &gmap);
                                let total2: f32 = get_sample_pair_total_ibd(blk2, &gmap);
                                total1_vec.push(total1);
                                total2_vec.push(total2);
                            }
                        }
                    }
                }
                false => {
                    for ibd in [&mut ibd1, &mut ibd2] {
                        if !ibd.is_sorted_by_haplotypes() {
                            ibd.sort_by_haplotypes();
                        }
                        ibd.infer_ploidy();
                    }
                    let it1 = ibd1
                        .as_slice()
                        .linear_group_by_key(|seg| seg.haplotype_pair());
                    let it2 = ibd2
                        .as_slice()
                        .linear_group_by_key(|seg| seg.haplotype_pair());
                    for e in it1.merge_join_by(it2, |blk1, blk2| {
                        blk1[0].haplotype_pair().cmp(&blk2[0].haplotype_pair())
                    }) {
                        // haplotype pairs, segments are not overlapping
                        // only need sum them up per haplotype pair
                        match e {
                            itertools::EitherOrBoth::Left(blk1) => {
                                let total1: f32 =
                                    blk1.iter().map(|seg| seg.get_seg_len_cm(&gmap)).sum();
                                total1_vec.push(total1);
                                total2_vec.push(0.0f32);
                            }
                            itertools::EitherOrBoth::Right(blk2) => {
                                let total2: f32 =
                                    blk2.iter().map(|seg| seg.get_seg_len_cm(&gmap)).sum();
                                total2_vec.push(total2);
                                total1_vec.push(0.0f32);
                            }
                            itertools::EitherOrBoth::Both(blk1, blk2) => {
                                let total1: f32 =
                                    blk1.iter().map(|seg| seg.get_seg_len_cm(&gmap)).sum();
                                total1_vec.push(total1);
                                let total2: f32 =
                                    blk2.iter().map(|seg| seg.get_seg_len_cm(&gmap)).sum();
                                total2_vec.push(total2);
                            }
                        }
                    }
                }
            }

            write_pair_total(
                total1_vec,
                total2_vec,
                "PairTotIbdA",
                "PairTotIbdB",
                out.with_extension("pairtotibdpq"),
            )?;
        }

        {
            // 3. population-level total ibd
            let binwidth = 0.05f32;
            let mut bincenters = vec![];
            let mut totals1 = vec![];
            let mut totals2 = vec![];

            for seg in ibd1.as_slice() {
                let cm = seg.get_seg_len_cm(&gmap);
                let idx = (cm / binwidth) as usize;
                if idx >= totals1.len() {
                    totals1.resize(idx + 1, 0.0f32);
                }
                totals1[idx] += cm;
            }

            for seg in ibd2.as_slice() {
                let cm = seg.get_seg_len_cm(&gmap);
                let idx = (cm / binwidth) as usize;
                if idx >= totals2.len() {
                    totals2.resize(idx + 1, 0.0f32);
                }
                totals2[idx] += cm;
            }

            // keep two vector same size
            let mut sz = totals1.len();
            if sz < totals2.len() {
                sz = totals2.len();
            }

            totals1.resize(sz, 0.0);
            totals2.resize(sz, 0.0);

            // make bin center vector
            for idx in 1..sz {
                let c = (idx as f32 + 0.5) * binwidth;
                bincenters.push(c);
            }

            // write csv file
            let mut f = std::fs::File::create(out.with_extension("poptotibdcsv"))?;
            use std::io::Write;
            writeln!(f, "BinCenter,PopTotIbdA,PopTotIbdB")?;
            for ((bc, tot1), tot2) in bincenters.iter().zip(totals1.iter()).zip(totals2.iter()) {
                writeln!(f, "{bc},{tot1},{tot2}")?;
            }
        }
    }
    Ok(())
}

/// calculate sample-pair total IBD in cM
///
/// Assumes (1) all segments are from the same sample pair;
/// (2) they are sorted by coordinates
///
/// When there are overlapping segments, merged IBD segment length is
/// calculated and accumulated.
fn get_sample_pair_total_ibd(blk: &[IbdSeg], gmap: &GeneticMap) -> f32 {
    let mut total = 0.0f32;
    let mut cur_seg = blk[0];
    blk.iter().skip(1).for_each(|seg| {
        if seg.s > cur_seg.e {
            total += cur_seg.get_seg_len_cm(gmap);
            cur_seg.s = seg.s;
            cur_seg.e = seg.e;
        } else {
            cur_seg.e = seg.e;
        }
    });
    if cur_seg.s < cur_seg.e {
        total += cur_seg.get_seg_len_cm(gmap);
    }

    total
}
