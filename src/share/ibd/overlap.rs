use super::ibdset::*;
use crate::container::intervals::Intervals;
use crate::container::intervaltree::IntervalTree;
use crate::genome::GenomeInfo;
use crate::gmap::GeneticMap;
use std::io::BufWriter;
use std::path::PathBuf;

pub struct IbdOverlapAnalyzer<'a> {
    ibd1: &'a IbdSet<'a>,
    ibd2: &'a IbdSet<'a>,
    prefix_for_details: Option<&'a PathBuf>,
    ignore_hap: bool,
}

impl<'a> IbdOverlapAnalyzer<'a> {
    pub fn new(
        ibd1: &'a IbdSet<'a>,
        ibd2: &'a IbdSet<'a>,
        ignore_hap: bool,
        prefix_for_details: Option<&'a PathBuf>,
    ) -> Self {
        if ignore_hap {
            assert!(ibd1.is_sorted_by_samples());
            assert!(ibd2.is_sorted_by_samples());
        } else {
            assert!(ibd1.is_sorted_by_haplotypes());
            assert!(ibd2.is_sorted_by_haplotypes());
        }
        assert!(ibd1.has_same_individuals(ibd2));
        Self {
            ibd1,
            ibd2,
            prefix_for_details,
            ignore_hap,
        }
    }

    fn calc_overlap_rates<'b>(
        &self,
        ibd1: &'b IbdSet<'a>,
        ibd2: &'b IbdSet<'a>,
        len_ranges: Option<&[f32]>,
        swap12: bool,
    ) -> (Vec<f64>, Vec<f64>, f64, f64) {
        let gmap = self.ibd1.get_gmap();
        let ginfo = self.ibd1.get_ginfo();
        let sample_names = self.ibd1.get_inds().v();
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        let mut counters = vec![0usize; len_ranges.len()];

        // define variance in terms of moments expression
        // mean = sum(x)/n
        // sd = sqrt(var)
        // var = sum(x^2)/n - (sum(x)/n)^
        let mut ratio_sums: Vec<f64> = vec![0.0f64; len_ranges.len()];
        let mut ratio_square_sums: Vec<f64> = vec![0.0f64; len_ranges.len()];
        let mut gw_counters = 0usize;
        let mut gw_ratio_sum = 0.0f64;
        let mut gw_ratio_square_sum = 0.0f64;
        let mut tree = IntervalTree::<u32, ()>::new(100);
        let mut itvs = Intervals::new();
        let mut detail_file = match (self.prefix_for_details, swap12) {
            (None, _) => None,
            (Some(prefix_for_detials), true) => Some(
                std::fs::File::create(prefix_for_detials.with_extension("abyb"))
                    .map(|f| BufWriter::with_capacity(100000000, f))
                    .unwrap(),
            ),
            (Some(prefix_for_detials), false) => Some(
                std::fs::File::create(prefix_for_detials.with_extension("bbya"))
                    .map(|f| BufWriter::with_capacity(100000000, f))
                    .unwrap(),
            ),
        };

        let blkpair_iter = if swap12 {
            IbdSetBlockPairIter::new(ibd2, ibd1, self.ignore_hap)
        } else {
            IbdSetBlockPairIter::new(ibd1, ibd2, self.ignore_hap)
        };

        for (blk1, blk2) in blkpair_iter {
            let mut gw_total_a = 0.0;
            let mut gw_total_intersect = 0.0;
            let pair_info = match (blk1, blk2) {
                (Some(blk1), Some(blk2)) => {
                    // for output details
                    let (sample1, hapid1, sample2, hapid2) = {
                        let (i, j, m, n) = blk1[0].haplotype_pair();
                        let s1 = &sample_names[i as usize];
                        let s2 = &sample_names[j as usize];
                        let h1 = match m {
                            0 => 1,
                            1 => 2,
                            _ => 0,
                        };
                        let h2 = match n {
                            0 => 1,
                            1 => 2,
                            _ => 0,
                        };
                        (s1, h1, s2, h2)
                    };
                    if self.ignore_hap {
                        // when ignore hap, we need to flattend segments
                        itvs.clear();
                        itvs.extend_from_iter(blk2.iter().map(|x| (x.s, x.e)));
                        itvs.merge();
                        let it = itvs.iter().map(|x| (x.clone(), ()));
                        tree.clear_and_fill_with_iter(it);
                    } else {
                        let it = blk2.iter().map(|seg| (seg.s..seg.e, ()));
                        tree.clear_and_fill_with_iter(it);
                    };

                    // flatten a
                    itvs.clear();
                    itvs.extend_from_iter(blk1.iter().map(|x| (x.s, x.e)));
                    if self.ignore_hap {
                        itvs.merge();
                    } else {
                        itvs.sort();
                    }

                    for a in itvs.iter() {
                        let s1 = a.start;
                        let e1 = a.end;
                        let cm = gmap.get_cm_len(s1, e1);

                        gw_total_a += cm;

                        // total intersected length of this a segment by any b segments
                        let mut intersect = 0.0f32;
                        for b in tree.query(a.clone()) {
                            let b = &b.range;
                            let s2 = b.start;
                            let e2 = b.end;
                            let mut s = s1;
                            if s2 > s {
                                s = s2;
                            }
                            let mut e = e1;
                            if e2 < e {
                                e = e2;
                            }

                            intersect += gmap.get_cm_len(s, e);
                        }
                        gw_total_intersect += intersect;

                        if cm < len_ranges[0] {
                            // if two short no need to compare
                            continue;
                        }
                        // length category
                        let which = len_ranges.partition_point(|x| *x <= cm) - 1;
                        // update the summary vectors
                        counters[which] += 1;
                        let ratio = (intersect / cm) as f64;
                        ratio_sums[which] += ratio;
                        ratio_square_sums[which] += ratio * ratio;

                        // Write a record to the detail file for each subject segment
                        // s1,h1,s2,h2,chr,astart, aend, overlap_ratio, interval_lwr
                        if let Some(detail_file) = detail_file.as_mut() {
                            use std::io::Write;
                            let interval_lwr = len_ranges[which];
                            let (_ichr, chrname, astart) = ginfo.to_chr_pos(a.start);
                            let aend = astart + (a.end - a.start);
                            write!(detail_file, "{sample1}\t{hapid1}\t{sample2}\t{hapid2}\t{chrname}\t{astart}\t{aend}\t{interval_lwr}\t{:.4}\n", intersect/cm ).unwrap();
                        }
                    }
                    Some((sample1, hapid1, sample2, hapid2))
                }
                (None, Some(_blk2)) => {
                    // ignore
                    None
                }
                (Some(blk1), None) => {
                    // for output details
                    let (sample1, hapid1, sample2, hapid2) = {
                        let (i, j, m, n) = blk1[0].haplotype_pair();
                        let s1 = &sample_names[i as usize];
                        let s2 = &sample_names[j as usize];
                        let h1 = match m {
                            0 => 1,
                            1 => 2,
                            _ => 0,
                        };
                        let h2 = match n {
                            0 => 1,
                            1 => 2,
                            _ => 0,
                        };
                        (s1, h1, s2, h2)
                    };
                    // flatten a
                    itvs.clear();
                    itvs.extend_from_iter(blk1.iter().map(|x| (x.s, x.e)));
                    if self.ignore_hap {
                        itvs.merge();
                    } else {
                        itvs.sort();
                    }
                    for a in itvs.iter() {
                        let s1 = a.start;
                        let e1 = a.end;
                        let cm = gmap.get_cm_len(s1, e1);

                        gw_total_a += cm;

                        if cm < len_ranges[0] {
                            // if two short no need to compare
                            continue;
                        }
                        // length category
                        let which = len_ranges.partition_point(|x| *x <= cm) - 1;
                        // total intersected length of this a segment by any b segments
                        let intersect = 0.0f32;
                        // update the summary vectors
                        counters[which] += 1;

                        let ratio = (intersect / cm) as f64;
                        ratio_sums[which] += ratio;
                        ratio_square_sums[which] += ratio * ratio;
                        // Write a record to the detail file for each subject segment
                        // s1,h1,s2,h2,chr,astart, aend, overlap_ratio, interval_lwr
                        if let Some(detail_file) = detail_file.as_mut() {
                            use std::io::Write;
                            let interval_lwr = len_ranges[which];
                            let (_ichr, chrname, astart) = ginfo.to_chr_pos(a.start);
                            let aend = astart + (a.end - a.start);
                            write!(detail_file, "{sample1}\t{hapid1}\t{sample2}\t{hapid2}\t{chrname}\t{astart}\t{aend}\t{interval_lwr}\t{:.4}\n", intersect/cm ).unwrap();
                        }
                    }
                    Some((sample1, hapid1, sample2, hapid2))
                }
                _ => None,
            };
            if let Some((sample1, hapid1, sample2, hapid2)) = pair_info {
                gw_counters += 1;
                let gw_ratio = (gw_total_intersect / gw_total_a) as f64;
                gw_ratio_sum += gw_ratio;
                gw_ratio_square_sum += gw_ratio * gw_ratio;
                // Write a record to the detail file for a sample-pair with non-zero IBD sharig in
                // ibd1
                if let Some(detail_file) = detail_file.as_mut() {
                    use std::io::Write;
                    write!(
                    detail_file,
                    "{sample1}\t{hapid1}\t{sample2}\t{hapid2}\tgenome_wide\t-1\t-1\t-1\t{:.4}\n",
                    gw_total_intersect / gw_total_a
                )
                    .unwrap();
                }
            }
        }

        // average out
        //
        // ratio means
        for (n, tot) in counters.iter().zip(ratio_sums.iter_mut()) {
            *tot /= *n as f64;
        }
        let ratio_means = ratio_sums; // average above        let ratio_square_means = ratio_square_sums; // average above

        // ratio stds
        for ((n, mean), tot) in counters
            .iter()
            .zip(ratio_means.iter())
            .zip(ratio_square_sums.iter_mut())
        {
            *tot /= *n as f64; // sq mean
            *tot -= mean * mean; // var
            *tot = tot.sqrt(); // std
        }
        let ratio_square_stds = ratio_square_sums;

        // gw mean
        let gw_ratio_mean = gw_ratio_sum / gw_counters as f64;
        // gw std
        gw_ratio_square_sum /= gw_counters as f64; // sq mean
        gw_ratio_square_sum -= gw_ratio_mean * gw_ratio_mean; // var
        let gw_ratio_std = gw_ratio_square_sum.sqrt();

        (ratio_means, ratio_square_stds, gw_ratio_mean, gw_ratio_std)
    }

    pub fn analzyze(
        &self,
        len_ranges: Option<&[f32]>,
        window_gw_bp: Option<(u32, u32)>,
    ) -> IbdOverlapResult {
        let mut ibd1 = self.ibd1;
        let mut ibd2 = self.ibd2;

        // empty
        let mut ibd1_win = IbdSet::new(
            self.ibd1.get_gmap(),
            self.ibd1.get_ginfo(),
            self.ibd1.get_inds(),
        );
        let mut ibd2_win = IbdSet::new(
            self.ibd2.get_gmap(),
            self.ibd2.get_ginfo(),
            self.ibd2.get_inds(),
        );

        if let Some((win_start, win_end)) = window_gw_bp {
            for seg_orig in self.ibd1.as_slice() {
                let mut seg = *seg_orig;
                if win_start > seg.s {
                    seg.s = win_start;
                }
                if win_end < seg.e {
                    seg.e = win_end;
                }
                if seg.s < seg.e {
                    ibd1_win.add(seg);
                }
            }
            for seg_orig in self.ibd2.as_slice() {
                let mut seg = *seg_orig;
                if win_start > seg.s {
                    seg.s = win_start;
                }
                if win_end < seg.e {
                    seg.e = win_end;
                }
                if seg.s < seg.e {
                    ibd2_win.add(seg);
                }
            }
            ibd1 = &ibd1_win;
            ibd2 = &ibd2_win;
        }

        let swap = false;
        let (a_by_b_means, a_by_b_stds, a_by_b_gw_mean, a_by_b_gw_std) =
            self.calc_overlap_rates(ibd1, ibd2, len_ranges, swap);
        let swap = true;
        let (b_by_a_means, b_by_a_stds, b_by_a_gw_mean, b_by_a_gw_std) =
            self.calc_overlap_rates(ibd1, ibd2, len_ranges, swap);
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        IbdOverlapResult {
            bins: len_ranges.into(),
            a_by_b_means,
            a_by_b_stds,
            a_by_b_gw_mean,
            a_by_b_gw_std,
            b_by_a_means,
            b_by_a_stds,
            b_by_a_gw_mean,
            b_by_a_gw_std,
        }
    }

    pub fn get_paired_totalibd(self) {}
}

pub struct IbdOverlapResult {
    bins: Vec<f32>,
    a_by_b_means: Vec<f64>,
    a_by_b_stds: Vec<f64>,
    a_by_b_gw_mean: f64,
    a_by_b_gw_std: f64,
    b_by_a_means: Vec<f64>,
    b_by_a_stds: Vec<f64>,
    b_by_a_gw_mean: f64,
    b_by_a_gw_std: f64,
}

impl IbdOverlapResult {
    pub fn to_csv(&self, p: impl AsRef<std::path::Path>) {
        use std::io::Write;
        let mut file = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .unwrap();

        write!(file, "LenBinStart,RateAOverlapByBMean,RateAOverlapByBStd,RateBOverlapByAMean,RateBOverlapByAStd\n").unwrap();
        for ((((bin, ab_mean), ab_std), ba_mean), ba_std) in self
            .bins
            .iter()
            .zip(self.a_by_b_means.iter())
            .zip(self.a_by_b_stds.iter())
            .zip(self.b_by_a_means.iter())
            .zip(self.b_by_a_stds.iter())
        {
            write!(file, "{bin},{ab_mean},{ab_std},{ba_mean},{ba_std}\n").unwrap();
        }
        write!(
            file,
            "genome_wide,{},{},{},{}\n",
            self.a_by_b_gw_mean, self.a_by_b_gw_std, self.b_by_a_gw_mean, self.b_by_a_gw_std
        )
        .unwrap();
    }
}

pub fn write_per_winddow_overlap_res(
    windows: &[(u32, u32)],
    ov_res_slice: &[IbdOverlapResult],
    ginfo: &GenomeInfo,
    gmap: &GeneticMap,
    p: impl AsRef<std::path::Path>,
) {
    use std::io::Write;
    let mut file = std::fs::File::create(p.as_ref())
        .map(BufWriter::new)
        .unwrap();

    // write csv header
    write!(file, "Chrom,StartBp,EndBp,StartCm,EndCm,GwStartBp,GwEndBp,GwStartCm,GwEndCm,LenBinStart,AOvByB,BOvByA\n").unwrap();
    for (&(win_start, win_end), ov_res) in windows.iter().zip(ov_res_slice.iter()) {
        let (idx_s, chr_s, chr_pos_s) = ginfo.to_chr_pos(win_start);
        let gwcm_s = gmap.get_cm(win_start);
        let chrcm_s = gwcm_s - gmap.get_cm(ginfo.gwstarts[idx_s]);
        let (idx_e, _chr_e, chr_pos_e) = ginfo.to_chr_pos(win_end);
        let gwcm_e = gmap.get_cm(win_end - 1);
        let chrcm_e = gwcm_e - gmap.get_cm(ginfo.gwstarts[idx_e]);
        assert_eq!(idx_s, idx_e);

        for ((bin, ab), ba) in ov_res
            .bins
            .iter()
            .zip(ov_res.a_by_b_means.iter())
            .zip(ov_res.b_by_a_means.iter())
        {
            write!(file, "{chr_s},{chr_pos_s},{chr_pos_e},{chrcm_s},{chrcm_e},{win_start},{win_end},{gwcm_s},{gwcm_e},{bin},{ab},{ba}\n").unwrap();
        }
        write!(
            file,
            "{chr_s},{chr_pos_s},{chr_pos_e},{chrcm_s},{chrcm_e},{win_start},{win_end},{gwcm_s},{gwcm_e},genome_wide,{},{}\n",
            ov_res.a_by_b_gw_mean, ov_res.b_by_a_gw_mean
        )
        .unwrap();
    }
}
