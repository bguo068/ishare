use super::ibdset::*;
use crate::container::intervals::Intervals;
use crate::container::intervaltree::IntervalTree;
use std::io::BufWriter;

pub struct IbdOverlapAnalyzer<'a> {
    ibd1: &'a IbdSet<'a>,
    ibd2: &'a IbdSet<'a>,
    ignore_hap: bool,
}

impl<'a> IbdOverlapAnalyzer<'a> {
    pub fn new(ibd1: &'a IbdSet<'a>, ibd2: &'a IbdSet<'a>, ignore_hap: bool) -> Self {
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
            ignore_hap,
        }
    }

    fn calc_overlap_rates(&self, len_ranges: Option<&[f32]>, swap12: bool) -> Vec<f32> {
        let gmap = self.ibd1.get_gmap();
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        let mut counters = vec![0usize; len_ranges.len()];
        let mut ratio_sums: Vec<f32> = vec![0.0f32; len_ranges.len()];
        let mut tree = IntervalTree::<u32, ()>::new(100);
        let mut itvs = Intervals::new();

        let blkpair_iter = if swap12 {
            IbdSetBlockPairIter::new(&self.ibd2, &self.ibd1, self.ignore_hap)
        } else {
            IbdSetBlockPairIter::new(&self.ibd1, &self.ibd2, self.ignore_hap)
        };

        for (blk1, blk2) in blkpair_iter {
            match (blk1, blk2) {
                (Some(blk1), Some(blk2)) => {
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
                        if cm < len_ranges[0] {
                            // if two short no need to compare
                            continue;
                        }
                        // length category
                        let which = len_ranges.partition_point(|x| *x <= cm) - 1;
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
                        // update the summary vectors
                        counters[which] += 1;
                        ratio_sums[which] += intersect / cm;
                    }
                }
                (None, Some(_blk2)) => {
                    // ignore
                }
                (Some(blk1), None) => {
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
                        ratio_sums[which] = intersect / cm;
                    }
                }
                _ => {}
            }
        }

        // average out
        for (n, tot) in counters.iter().zip(ratio_sums.iter_mut()) {
            *tot /= *n as f32;
        }

        ratio_sums
    }

    pub fn analzyze(&self, len_ranges: Option<&[f32]>) -> IbdOverlapResult {
        let a_by_b = self.calc_overlap_rates(len_ranges, false);
        let b_by_a = self.calc_overlap_rates(len_ranges, true);
        let len_ranges = match len_ranges {
            Some(x) => x,
            None => &[0.0f32],
        };
        IbdOverlapResult {
            bins: len_ranges.into(),
            a_by_b,
            b_by_a,
        }
    }

    pub fn get_paired_totalibd(self) {}
}

pub struct IbdOverlapResult {
    bins: Vec<f32>,
    a_by_b: Vec<f32>,
    b_by_a: Vec<f32>,
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
    }
}
