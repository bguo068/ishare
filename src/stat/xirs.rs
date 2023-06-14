use crate::io::*;
use crate::share::ibd::{ibdseg::IbdSeg, ibdset::*};

use statrs::distribution::{ChiSquared, ContinuousCDF};
/// Calcuate Xirs stats for each SNP
///
/// Input:
///
/// A matrix of binary IBD status with rows corresponding to SNPs and columns
/// corresponding to isolate pairs
///
/// Step:
/// 1. Subtract the column mean from all rows to account for the amount of
/// relatedness between each pair.
/// 2. Subtract the row mean from each row and divide by the square root of
/// pi(1-pi), where pi is the population allele frequency of SNP i. This adjusts
/// for differences in SNP allele frequencies, which can affect the ability to
/// detect IBD
/// 3. Calculate row sums and divide these values by the square root of the
/// number of pairs (raw summary statistics)
/// 4. Normalize the raw summary statistics genome-wisely by binning all SNPs
/// into 100 equally sized bins partitioned on allele frequencies and then we
/// subtracted the mean and divided by the standard deviation of all values
/// within each bin (z-score)
/// 5. Square z-score (as negative is hard to interpret) to convert to
/// a new stats (Xir,s) that follows a chi-squared distribution with 1 degree of freedom.
///
///
/// Reference:
///
/// Henden et al 2018 [https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007279]
pub struct XirsBuilder<'a> {
    afrq: Vec<f64>,
    site_pos: Vec<u32>,
    ibd: &'a IbdSet<'a>,
    nhap: usize,
    cm: Vec<f64>,   // intermediate variables (columns means in step 1)
    rm: Vec<f64>,   // intermediate variables (row means, in step 2)
    rs: Vec<f64>,   // intermediate variables (row sum, in step 3)
    raw: Vec<f64>,  // un-normalized raw stats
    norm: Vec<f64>, // normalized raw stats
    xirs: Vec<f64>,
    pval: Vec<f64>,
}

impl<'a> XirsBuilder<'a> {
    pub fn new(afrq: Vec<f64>, site_pos: Vec<u32>, ibd: &'a IbdSet<'a>) -> Self {
        let mut site_pos = site_pos;
        // sort and check positions order and range
        site_pos.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!(
            site_pos
                .iter()
                .zip(site_pos.iter().skip(1))
                .all(|(x, y)| x.partial_cmp(y).unwrap().is_lt()),
            "unique positions"
        );
        assert!(*site_pos.last().unwrap() <= ibd.get_ginfo().get_total_len_bp());

        // check p and calcuate pq_sqrt
        assert!(afrq.iter().all(|x| (*x < 1.0) && (*x > 0.0)));

        assert_eq!(site_pos.len(), afrq.len());

        assert!(ibd.is_sorted_by_haplotypes());
        assert!(ibd.iter().all(|x| x.is_valid() && (!x.is_from_merge())));

        let nhap = ibd.get_inds().v().len() * 2;
        let nsites = site_pos.len();
        let npairs = nhap * (nhap - 1) / 2;
        Self {
            afrq,
            site_pos,
            ibd,
            nhap,
            cm: Vec::with_capacity(npairs),
            rm: Vec::with_capacity(nsites),
            rs: Vec::with_capacity(nsites),
            raw: Vec::with_capacity(nsites),
            norm: Vec::with_capacity(nsites),
            xirs: Vec::with_capacity(nsites),
            pval: Vec::with_capacity(nsites),
        }
    }

    pub fn get_total_num_pairs(&self) -> usize {
        self.nhap * (self.nhap - 1) / 2
    }
    pub fn get_total_num_sites(&self) -> u32 {
        self.site_pos.len() as u32
    }

    // /// array idx to upper matrix idx
    // fn x2ij(x: usize) -> (u32, u32) {
    //     let i = ((2.0 * x as f64 + 0.25).sqrt() + 0.5) as usize;
    //     // j here need to be wide enough to avoid overflow during multiplification
    //     let j = x - i * (i - 1) / 2;
    //     (i as u32, j as u32)
    // }
    fn ij2x(i: u32, j: u32) -> usize {
        let j = j as usize;
        let i = i as usize;
        (i - 1) * i / 2 + j
    }

    fn seg_2_x(seg: &IbdSeg) -> usize {
        let (id1, hap1, id2, hap2) = seg.haplotype_pair();
        let hap_i = (id1 << 1) + (hap1 as u32);
        let hap_j = (id2 << 1) + (hap2 as u32);
        let col = Self::ij2x(hap_i, hap_j);
        col
    }

    fn blk_2_x(blk: &[IbdSeg]) -> usize {
        let seg = &blk[0];
        Self::seg_2_x(seg)
    }

    fn calculate_cm(&mut self) {
        let p = self.site_pos.as_slice();
        let m = self.get_total_num_sites() as usize;
        let n = self.get_total_num_pairs();

        let mut sum = vec![0f64; n];

        for seg in self.ibd.iter() {
            let start = p.partition_point(|x| *x < seg.s);
            let end = start + p[start..].partition_point(|x| *x < seg.e);
            let col = Self::seg_2_x(seg);
            for _row in start..end {
                sum[col] += 1.0;
            }
        }

        self.cm.extend(sum.into_iter().map(|x| x / m as f64));
    }

    fn calculate_rm(&mut self) {
        let p = self.site_pos.as_slice();
        let m = self.get_total_num_sites() as usize;
        let n = self.get_total_num_pairs();

        let mut sum = vec![0f64; m];
        // for hap pairs that share at least one IBD segment
        for blk in IbdSetBlockIter::new(&self.ibd, false) {
            let col = Self::blk_2_x(blk);

            let mut non_ibd_start = 0usize;
            let mut non_ibd_end: usize;
            let mut start = 0usize;
            let mut end: usize;
            for seg in blk {
                start = start + p[start..].partition_point(|x| *x < seg.s);
                end = start + p[start..].partition_point(|x| *x < seg.e);
                // non-ibd region
                non_ibd_end = start;
                for row in non_ibd_start..non_ibd_end {
                    sum[row] += 0.0 - self.cm[col];
                }
                // ibd region
                for row in start..end {
                    sum[row] += 1.0 - self.cm[col];
                }
                non_ibd_start = end;
                start = end;
            }
            // trailing non-ibd region
            for row in non_ibd_start..p.len() {
                sum[row] += 0.0 - self.cm[col];
            }
        }
        // for hap pair that does not share IBD segments
        let it_all = (0..n).into_iter();
        let it_shared = IbdSetBlockIter::new(&self.ibd, false).map(Self::blk_2_x);
        use itertools::EitherOrBoth::*;
        let merged = itertools::merge_join_by(it_all, it_shared, |a, b| a.cmp(b));
        merged.for_each(|x| match x {
            Left(col) => {
                for row in 0..p.len() {
                    sum[row] += 0.0 - self.cm[col];
                }
            }
            _ => {}
        });

        self.rm.extend(sum.into_iter().map(|x| x / n as f64));
    }

    fn calculate_rs(&mut self) {
        let p = self.site_pos.as_slice();
        let m = self.get_total_num_sites() as usize;
        let n = self.get_total_num_pairs();
        let pqsqrt: Vec<f64> = self.afrq.iter().map(|p| (*p * (1.0 - *p)).sqrt()).collect();

        let mut sum = vec![0f64; m];
        // for hap pairs that share at least one IBD segment
        for blk in IbdSetBlockIter::new(&self.ibd, false) {
            let col = Self::blk_2_x(blk);

            let mut non_ibd_start = 0usize;
            let mut non_ibd_end: usize;
            let mut start = 0usize;
            let mut end: usize;
            for seg in blk {
                start = start + p[start..].partition_point(|x| *x < seg.s);
                end = start + p[start..].partition_point(|x| *x < seg.e);
                // non-ibd region
                non_ibd_end = start;
                for row in non_ibd_start..non_ibd_end {
                    sum[row] += (0.0 - self.cm[col] - self.rm[row]) / pqsqrt[row];
                }
                // ibd region
                for row in start..end {
                    sum[row] += (1.0 - self.cm[col] - self.rm[row]) / pqsqrt[row];
                }
                non_ibd_start = end;
                start = end;
            }
            // trailing non-ibd region
            for row in non_ibd_start..p.len() {
                sum[row] += (0.0 - self.cm[col] - self.rm[row]) / pqsqrt[row];
            }
        }
        // for hap pair that does not share IBD segments
        let it_all = (0..n).into_iter();
        let it_shared = IbdSetBlockIter::new(&self.ibd, false).map(Self::blk_2_x);
        use itertools::EitherOrBoth::*;
        let merged = itertools::merge_join_by(it_all, it_shared, |a, b| a.cmp(b));
        merged.for_each(|x| match x {
            Left(col) => {
                for row in 0..p.len() {
                    sum[row] += (0.0 - self.cm[col] - self.rm[row]) / pqsqrt[row];
                }
            }
            _ => {}
        });
        self.rs = sum;
    }

    fn calculate_raw(&mut self) {
        if self.rs.len() != self.get_total_num_sites() as usize {
            self.calculate_rs();
        }
        for rs in self.rs.iter() {
            self.raw
                .push(*rs / (self.get_total_num_pairs() as f64).sqrt());
        }
    }

    // 4. Normalize the raw summary statistics genome-wisely by binning all SNPs
    // into 100 equally sized bins partitioned on allele frequencies and then we
    // subtracted the mean and divided by the standard deviation of all values
    // within each bin (z-score)
    fn calculate_norm(&mut self) {
        if self.raw.len() != self.get_total_num_sites() as usize {
            self.calculate_raw();
        }

        // vec of 3-tuple (orig_idx, freq, class)
        let mut annot_freq: Vec<(usize, f64, usize)> = self
            .afrq
            .iter()
            .enumerate()
            .map(|(i, p)| (i, *p, 0))
            .collect();

        // sort by allele frequencies
        annot_freq.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // assign frquency bin class
        // the first `extra` bins have `binsz + 1` members
        // the rest `100 - extra` bins have `binsz` members.
        let n = self.get_total_num_sites() as usize;
        let binsz = n / 100;
        let extra = n % 100;
        let mut s = 0usize;
        let mut class = 0;
        while s < n {
            let e = s + binsz + (class < extra) as usize;
            annot_freq[s..e].iter_mut().for_each(|x| x.2 = class);
            class += 1;
            s = e;
        }

        // revert back to the original order
        annot_freq.sort_by_key(|x| x.0);

        // summary per class
        let mut run_sums = Vec::<f64>::with_capacity(100);
        let mut run_cnts = Vec::<usize>::with_capacity(100);

        run_sums.resize(100, 0.0);
        run_cnts.resize(100, 0);
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw.iter()) {
            run_cnts[*which] += 1;
            run_sums[*which] += *raw_val;
        }

        let means: Vec<f64> = run_sums
            .iter()
            .zip(run_cnts.iter())
            .map(|(s, n)| *s / (*n) as f64)
            .collect();

        run_sums.resize(100, 0.0);
        run_cnts.resize(100, 0);
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw.iter()) {
            run_cnts[*which] += 1;
            let diff = *raw_val - means[*which];
            run_sums[*which] += diff * diff;
        }

        let stds: Vec<f64> = run_sums
            .iter()
            .zip(run_cnts.iter())
            .map(|(s, n)| (*s / (*n) as f64).sqrt())
            .collect();

        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw.iter()) {
            let which = *which;
            let m = means[which];
            let std = stds[which];
            self.norm.push((*raw_val - m) / std);
        }
    }

    fn calculate_xirs(&mut self) {
        if self.norm.len() != self.get_total_num_sites() as usize {
            self.calculate_norm();
        }

        self.xirs.clear();
        // squre of normalized stats
        self.xirs.extend(self.norm.iter().map(|x| *x * *x));
    }

    fn calculate_pval(&mut self) {
        if self.xirs.len() != self.get_total_num_sites() as usize {
            self.calculate_xirs();
        }
        let chisq = ChiSquared::new(1.0).unwrap();

        self.pval.clear();
        // calculate pvalue using cdf function
        self.pval
            .extend(self.xirs.iter().map(|x| 1.0 - chisq.cdf(*x)));
    }

    pub fn finish(&mut self) -> XirsResult {
        println!("calculate 1. cm");
        self.calculate_cm();
        println!("calculate 2. rm");
        self.calculate_rm();
        println!("calculate 3. rs");
        self.calculate_rs();
        println!("calculate 4- others");
        self.calculate_raw();
        self.calculate_norm();
        self.calculate_xirs();
        self.calculate_pval();

        let n = self.get_total_num_sites() as usize;
        let mut chr_id = Vec::<u32>::with_capacity(n);
        let mut chr_pos = Vec::<u32>::with_capacity(n);
        let gw_pos = self.site_pos.clone();

        for gwpos in gw_pos.iter() {
            let (chid, _, chrpos) = self.ibd.get_ginfo().to_chr_pos(*gwpos);
            chr_id.push(chid as u32);
            chr_pos.push(chrpos);
        }
        let xirs = self.xirs.clone();
        let pval = self.pval.clone();

        XirsResult {
            chr_id,
            chr_pos,
            gw_pos,
            xirs,
            pval,
        }
    }
}

pub struct XirsResult {
    pub chr_id: Vec<u32>,
    pub chr_pos: Vec<u32>,
    pub gw_pos: Vec<u32>,
    pub xirs: Vec<f64>,
    pub pval: Vec<f64>,
}

impl IntoParquet for XirsResult {
    fn into_parquet(&mut self, p: impl AsRef<std::path::Path>) {
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::arrow_writer::ArrowWriter;
        use parquet::file::properties::WriterProperties;
        use std::fs::File;
        use std::mem::take;

        let batch = RecordBatch::try_from_iter(vec![
            ("ChrId", take(&mut self.chr_id).into_arrow_array()),
            ("ChrPos", take(&mut self.chr_pos).into_arrow_array()),
            ("GwPos", take(&mut self.gw_pos).into_arrow_array()),
            ("Xirs", take(&mut self.xirs).into_arrow_array()),
            ("Pval", take(&mut self.pval).into_arrow_array()),
        ])
        .unwrap();

        let file = File::create(p.as_ref()).unwrap();
        // Default writer properties
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();

        writer.write(&batch).expect("Writing batch");

        // writer must be closed to write footer
        writer.close().unwrap();
    }
}
