use super::ibdset::*;
use crate::container::intervals::Intervals;
use crate::container::intervaltree::IntervalTree;
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

    fn calc_overlap_rates(&self, len_ranges: Option<&[f32]>, swap12: bool) -> (Vec<f32>, f32) {
        let gmap = self.ibd1.get_gmap();
        let ginfo = self.ibd1.get_ginfo();
        let sample_names = self.ibd1.get_inds().v();
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        let mut counters = vec![0usize; len_ranges.len()];
        let mut ratio_sums: Vec<f32> = vec![0.0f32; len_ranges.len()];
        let mut gw_counters = 0usize;
        let mut gw_ratio_sums = 0.0f32;
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
            IbdSetBlockPairIter::new(&self.ibd2, &self.ibd1, self.ignore_hap)
        } else {
            IbdSetBlockPairIter::new(&self.ibd1, &self.ibd2, self.ignore_hap)
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
                        ratio_sums[which] += intersect / cm;

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

                        ratio_sums[which] += intersect / cm;
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
                gw_ratio_sums += gw_total_intersect / gw_total_a;
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
        for (n, tot) in counters.iter().zip(ratio_sums.iter_mut()) {
            *tot /= *n as f32;
        }

        (ratio_sums, gw_ratio_sums / gw_counters as f32)
    }

    pub fn analzyze(&self, len_ranges: Option<&[f32]>) -> IbdOverlapResult {
        let swap = false;
        let (a_by_b, a_by_b_gw) = self.calc_overlap_rates(len_ranges, swap);
        let swap = true;
        let (b_by_a, b_by_a_gw) = self.calc_overlap_rates(len_ranges, swap);
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        IbdOverlapResult {
            bins: len_ranges.into(),
            a_by_b,
            a_by_b_gw,
            b_by_a,
            b_by_a_gw,
        }
    }

    pub fn get_paired_totalibd(self) {}
}

pub struct IbdOverlapResult {
    bins: Vec<f32>,
    a_by_b: Vec<f32>,
    a_by_b_gw: f32,
    b_by_a: Vec<f32>,
    b_by_a_gw: f32,
}

impl IbdOverlapResult {
    pub fn to_csv(&self, p: impl AsRef<std::path::Path>) {
        use std::io::Write;
        let mut file = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .unwrap();

        write!(file, "LenBinStart,RateAOverlapByB,RateBOverlapByA\n").unwrap();
        for ((bin, ab), ba) in self
            .bins
            .iter()
            .zip(self.a_by_b.iter())
            .zip(self.b_by_a.iter())
        {
            write!(file, "{bin},{ab},{ba}\n").unwrap();
        }
        write!(file, "genome_wide,{},{}\n", self.a_by_b_gw, self.b_by_a_gw).unwrap();
    }
}
