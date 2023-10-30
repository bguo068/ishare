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
/// Henden et al 2018 <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007279>
pub struct XirsBuilder<'a> {
    afrq_i: Vec<f64>,
    pos_i: Vec<u32>,
    ibd: &'a IbdSet<'a>,
    npairs: usize,
    nsites: usize,
    xs_i: Vec<f64>,   // intermediate variables (columns means in step 1)
    cm_j: Vec<f64>,   // intermediate variables (columns means in step 1)
    rm_i: Vec<f64>,   // intermediate variables (row means, in step 2)
    rs_i: Vec<f64>,   // intermediate variables (row sum, in step 3)
    raw_i: Vec<f64>,  // un-normalized raw stats
    norm_i: Vec<f64>, // normalized raw stats
    xirs_i: Vec<f64>,
    pval_i: Vec<f64>,
}

impl<'a> XirsBuilder<'a> {
    pub fn new(afrq: Vec<f64>, site_pos: Vec<u32>, ibd: &'a IbdSet<'a>) -> Self {
        // sort and check positions order and range
        assert!(
            site_pos
                .iter()
                .zip(site_pos.iter().skip(1))
                .all(|(x, y)| { x.partial_cmp(y).unwrap().is_lt() }),
            "unique positions"
        );
        assert!(*site_pos.last().unwrap() <= ibd.get_ginfo().get_total_len_bp());

        // check p and calcuate pq_sqrt
        assert!(afrq.iter().all(|x| (*x < 1.0) && (*x > 0.0)));
        assert_eq!(site_pos.len(), afrq.len());
        assert!(ibd.is_sorted_by_haplotypes());
        assert!(ibd.iter().all(|x| x.is_valid() && (!x.is_from_merge())));

        let nsites = site_pos.len();
        let nhap = match ibd.get_ploidy_status() {
            IbdSetPloidyStatus::Haploid => ibd.get_inds().v().len(),
            IbdSetPloidyStatus::Diploid => ibd.get_inds().v().len() * 2,
            _ => {
                panic!("Xirs can only be calculated for two IbdSetPloidStatus: Haploid or Diploid");
            }
        };
        let npairs = nhap * (nhap - 1) / 2;
        Self {
            afrq_i: afrq,
            pos_i: site_pos,
            ibd,
            npairs,
            nsites,
            xs_i: Vec::with_capacity(nsites),
            cm_j: Vec::with_capacity(npairs),
            rm_i: Vec::with_capacity(nsites),
            rs_i: Vec::with_capacity(nsites),
            raw_i: Vec::with_capacity(nsites),
            norm_i: Vec::with_capacity(nsites),
            xirs_i: Vec::with_capacity(nsites),
            pval_i: Vec::with_capacity(nsites),
        }
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

    fn seg_2_x(&self, seg: &IbdSeg) -> usize {
        let (id1, hap1, id2, hap2) = seg.haplotype_pair();
        let (hap_i, hap_j) = match self.ibd.get_ploidy_status() {
            IbdSetPloidyStatus::Haploid => (id1, id2),
            IbdSetPloidyStatus::Diploid => ((id1 << 1) + (hap1 as u32), (id2 << 1) + (hap2 as u32)),
            _ => {
                panic!("Xirs can only be calculated for two IbdSetPloidStatus: Haploid or Diploid");
            }
        };

        let col = Self::ij2x(hap_i, hap_j);
        col
    }

    // fn blk_2_x(blk: &[IbdSeg]) -> usize {
    //     let seg = &blk[0];
    //     Self::seg_2_x(seg)
    // }

    fn calculate_cm_xs(&mut self) {
        let p = self.pos_i.as_slice();
        let m = self.nsites;
        let n = self.npairs;

        let mut xs = vec![0f64; m];
        let mut sum = vec![0f64; n];

        for seg in self.ibd.iter() {
            let start = p.partition_point(|x| *x < seg.s);
            let end = start + p[start..].partition_point(|x| *x < seg.e);
            let col = self.seg_2_x(seg);
            for row in start..end {
                sum[col] += 1.0;
                xs[row] += 1.0;
            }
        }

        self.cm_j.extend(sum.into_iter().map(|x| x / m as f64));
        self.xs_i = xs;
    }

    fn calculate_rm(&mut self) {
        // let p = self.pos_i.as_slice();
        // let m = self.nsites;
        let n = self.npairs;

        let cm_j_total: f64 = self.cm_j.iter().sum();

        let it = self.xs_i.iter().map(|x_i| (*x_i - cm_j_total) / (n as f64));
        self.rm_i.extend(it);
    }

    fn calculate_rs(&mut self) {
        // let p = self.pos_i.as_slice();
        // let m = self.nsites;
        let n = self.npairs;
        let pqsqrt: Vec<f64> = self
            .afrq_i
            .iter()
            .map(|p| (*p * (1.0 - *p)).sqrt())
            .collect();
        let cm_j_total: f64 = self.cm_j.iter().sum();

        let it = self
            .xs_i
            .iter()
            .zip(self.rm_i.iter())
            .zip(pqsqrt.iter())
            .map(|((xs_i, rm_i), pq_i)| {
                let mut val = *xs_i - cm_j_total - *rm_i * (n as f64);
                val /= pq_i;
                val
            });

        self.rs_i.extend(it);
    }

    fn calculate_raw(&mut self) {
        let m = self.nsites;
        if self.rs_i.len() != m {
            self.calculate_rs();
        }
        let n = self.npairs;
        for rs in self.rs_i.iter() {
            self.raw_i.push(*rs / (n as f64).sqrt());
        }
    }

    // 4. Normalize the raw summary statistics genome-wisely by binning all SNPs
    // into 100 equally sized bins partitioned on allele frequencies and then we
    // subtracted the mean and divided by the standard deviation of all values
    // within each bin (z-score)
    fn calculate_norm(&mut self) {
        let m = self.nsites;
        // let n = self.npairs;
        if self.raw_i.len() != m {
            self.calculate_raw();
        }

        // vec of 3-tuple (orig_idx, freq, class)
        let mut annot_freq: Vec<(usize, f64, usize)> = self
            .afrq_i
            .iter()
            .enumerate()
            .map(|(i, p)| (i, *p, 0))
            .collect();

        // sort by allele frequencies
        annot_freq.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // assign frquency bin class
        // the first `extra` bins have `binsz + 1` members
        // the rest `100 - extra` bins have `binsz` members.
        let binsz = m / 100;
        let extra = m % 100;
        let mut s = 0usize;
        let mut class = 0;
        while s < m {
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
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
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
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
            run_cnts[*which] += 1;
            let diff = *raw_val - means[*which];
            run_sums[*which] += diff * diff;
        }

        let stds: Vec<f64> = run_sums
            .iter()
            .zip(run_cnts.iter())
            .map(|(s, n)| (*s / (*n) as f64).sqrt())
            .collect();

        println!(
            "number of class with zero std: {}",
            stds.iter().map(|x| (*x < 1e-7) as u32).sum::<u32>()
        );

        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
            let which = *which;
            let m = means[which];
            let std = stds[which];
            self.norm_i.push((*raw_val - m) / std);
        }
    }

    fn calculate_xirs(&mut self) {
        let m = self.nsites;
        if self.norm_i.len() != m {
            self.calculate_norm();
        }

        self.xirs_i.clear();
        // squre of normalized stats
        self.xirs_i.extend(self.norm_i.iter().map(|x| *x * *x));
    }

    fn calculate_pval(&mut self) {
        let m = self.nsites;
        if self.xirs_i.len() != m {
            self.calculate_xirs();
        }
        let chisq = ChiSquared::new(1.0).unwrap();

        self.pval_i.clear();
        // calculate pvalue using cdf function
        self.pval_i
            .extend(self.xirs_i.iter().map(|x| 1.0 - chisq.cdf(*x)));
    }

    pub fn finish(&mut self) -> XirsResult {
        // println!("calculate 1. cm");
        self.calculate_cm_xs();
        // println!("calculate 2. rm");
        self.calculate_rm();
        // println!("calculate 3. rs");
        self.calculate_rs();
        // println!("calculate 4. raw");
        self.calculate_raw();
        // println!("calculate 5. norm");
        self.calculate_norm();
        // println!("calculate 6. xirs");
        self.calculate_xirs();
        // println!("calculate 7. xirs");
        self.calculate_pval();
        // println!("calculate 8. writing");

        // println!("cm: {}, {:?}..", self.cm_j.len(), &self.rm_i[1..5]);
        // println!("rm: {}, {:?}..", self.rm_i.len(), &self.rm_i[1..5]);
        // println!("rs: {}, {:?}..", self.rs_i.len(), &self.rs_i[1..5]);
        // println!("xs: {}, {:?}..", self.xs_i.len(), &self.xs_i[1..5]);
        // println!("raw: {}, {:?}..", self.raw_i.len(), &self.raw_i[1..5]);
        // println!("afreq: {}, {:?}..", self.afrq_i.len(), &self.afrq_i[1..5]);
        // println!("xris: {}, {:?}..", self.xirs_i.len(), &self.xirs_i[1..5]);

        let m = self.nsites;
        let mut chr_id = Vec::<u32>::with_capacity(m);
        let mut chr_pos = Vec::<u32>::with_capacity(m);
        let gw_pos = self.pos_i.clone();

        for gwpos in gw_pos.iter() {
            let (chid, _, chrpos) = self.ibd.get_ginfo().to_chr_pos(*gwpos);
            chr_id.push(chid as u32);
            chr_pos.push(chrpos);
        }
        let xirs = self.xirs_i.clone();
        let pval = self.pval_i.clone();
        let raw = self.raw_i.clone();

        XirsResult {
            chr_id,
            chr_pos,
            gw_pos,
            raw,
            xirs,
            pval,
        }
    }
}

pub struct XirsResult {
    pub chr_id: Vec<u32>,
    pub chr_pos: Vec<u32>,
    pub gw_pos: Vec<u32>,
    pub raw: Vec<f64>,
    pub xirs: Vec<f64>,
    pub pval: Vec<f64>,
}
pub struct XirsBuilder2<'a> {
    afrq_i: Vec<f64>,
    pos_i: Vec<u32>,
    ibd: &'a IbdSet<'a>,
    npairs: usize,
    nsites: usize,
    raw_i: Vec<f64>,  //
    norm_i: Vec<f64>, // normalized raw stats
    xirs_i: Vec<f64>,
    pval_i: Vec<f64>,
}

impl<'a> XirsBuilder2<'a> {
    pub fn new(afrq: Vec<f64>, site_pos: Vec<u32>, ibd: &'a IbdSet<'a>) -> Self {
        // sort and check positions order and range
        assert!(
            site_pos
                .iter()
                .zip(site_pos.iter().skip(1))
                .all(|(x, y)| { x.partial_cmp(y).unwrap().is_lt() }),
            "unique positions"
        );
        assert!(*site_pos.last().unwrap() <= ibd.get_ginfo().get_total_len_bp());

        // check p and calcuate pq_sqrt
        assert!(afrq.iter().all(|x| (*x < 1.0) && (*x > 0.0)));
        assert_eq!(site_pos.len(), afrq.len());
        assert!(ibd.is_sorted_by_haplotypes());
        assert!(ibd.iter().all(|x| x.is_valid() && (!x.is_from_merge())));

        let nsites = site_pos.len();
        let nhap = match ibd.get_ploidy_status() {
            IbdSetPloidyStatus::Haploid => ibd.get_inds().v().len(),
            IbdSetPloidyStatus::Diploid => ibd.get_inds().v().len() * 2,
            _ => {
                panic!("Xirs can only be calculated for two IbdSetPloidStatus: Haploid or Diploid");
            }
        };
        let npairs = nhap * (nhap - 1) / 2;
        Self {
            afrq_i: afrq,
            pos_i: site_pos,
            ibd,
            npairs,
            nsites,
            raw_i: Vec::with_capacity(nsites),
            norm_i: Vec::with_capacity(nsites),
            xirs_i: Vec::with_capacity(nsites),
            pval_i: Vec::with_capacity(nsites),
        }
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

    fn seg_2_x(&self, seg: &IbdSeg) -> usize {
        let (id1, hap1, id2, hap2) = seg.haplotype_pair();
        let (hap_i, hap_j) = match self.ibd.get_ploidy_status() {
            IbdSetPloidyStatus::Haploid => (id1, id2),
            IbdSetPloidyStatus::Diploid => ((id1 << 1) + (hap1 as u32), (id2 << 1) + (hap2 as u32)),
            _ => {
                panic!("Xirs can only be calculated for two IbdSetPloidStatus: Haploid or Diploid");
            }
        };

        let col = Self::ij2x(hap_i, hap_j);
        col
    }

    fn calculate_raw(&mut self) {
        let p = self.pos_i.as_slice();
        let m = self.nsites;
        let n = self.npairs;

        let mut c_j = vec![0f64; n];
        let mut r_i = vec![0f64; m];
        let mut x_rs_i = vec![0f64; m];
        let mut xx_rs_i = vec![0f64; m];
        let mut xc_rs_i = vec![0f64; m];

        // scan ibd to calcualte c_j first
        for seg in self.ibd.iter() {
            let start = p.partition_point(|x| *x < seg.s);
            let end = start + p[start..].partition_point(|x| *x < seg.e);
            let col = self.seg_2_x(seg);
            for _row in start..end {
                c_j[col] += 1.0;
            }
        }
        c_j.iter_mut().for_each(|x| *x /= m as f64);

        // calculate row aggregates
        for seg in self.ibd.iter() {
            let start = p.partition_point(|x| *x < seg.s);
            let end = start + p[start..].partition_point(|x| *x < seg.e);
            let col = self.seg_2_x(seg);
            for row in start..end {
                x_rs_i[row] += 1.0;
                xc_rs_i[row] += 1.0 * c_j[col];
                xx_rs_i[row] += 1.0;
            }
        }

        // calculate r_j
        let cj_sum: f64 = c_j.iter().sum();
        for i in 0..m {
            r_i[i] = (x_rs_i[i] - cj_sum) / n as f64;
        }

        for i in 0..m {
            let mut out = xx_rs_i[i];
            out += c_j.iter().map(|x| *x * *x).sum::<f64>();
            out += n as f64 * r_i[i] * r_i[i];
            out -= 2.0 * xc_rs_i[i];
            out -= 2.0 * r_i[i] * x_rs_i[i];
            out += 2.0 * r_i[i] * cj_sum;
            out /= self.afrq_i[i] * (1.0 - self.afrq_i[i]);
            self.raw_i.push(out);
        }
    }

    // 4. Normalize the raw summary statistics genome-wisely by binning all SNPs
    // into 100 equally sized bins partitioned on allele frequencies and then we
    // subtracted the mean and divided by the standard deviation of all values
    // within each bin (z-score)
    fn calculate_norm(&mut self) {
        let m = self.nsites;
        // let n = self.npairs;
        if self.raw_i.len() != m {
            self.calculate_raw();
        }

        // vec of 3-tuple (orig_idx, freq, class)
        let mut annot_freq: Vec<(usize, f64, usize)> = self
            .afrq_i
            .iter()
            .enumerate()
            .map(|(i, p)| (i, *p, 0))
            .collect();

        // sort by allele frequencies
        annot_freq.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // assign frquency bin class
        // the first `extra` bins have `binsz + 1` members
        // the rest `100 - extra` bins have `binsz` members.
        let binsz = m / 100;
        let extra = m % 100;
        let mut s = 0usize;
        let mut class = 0;
        while s < m {
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
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
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
        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
            run_cnts[*which] += 1;
            let diff = *raw_val - means[*which];
            run_sums[*which] += diff * diff;
        }

        let stds: Vec<f64> = run_sums
            .iter()
            .zip(run_cnts.iter())
            .map(|(s, n)| (*s / (*n) as f64).sqrt())
            .collect();

        println!(
            "number of class with zero std: {}",
            stds.iter().map(|x| (*x < 1e-7) as u32).sum::<u32>()
        );

        for ((_idx, _p, which), raw_val) in annot_freq.iter().zip(self.raw_i.iter()) {
            let which = *which;
            let m = means[which];
            let std = stds[which];
            self.norm_i.push((*raw_val - m) / std);
        }
    }

    fn calculate_xirs(&mut self) {
        let m = self.nsites;
        if self.norm_i.len() != m {
            self.calculate_norm();
        }

        self.xirs_i.clear();
        // squre of normalized stats
        self.xirs_i.extend(self.norm_i.iter().map(|x| *x * *x));
    }

    fn calculate_pval(&mut self) {
        let m = self.nsites;
        if self.xirs_i.len() != m {
            self.calculate_xirs();
        }
        let chisq = ChiSquared::new(1.0).unwrap();

        self.pval_i.clear();
        // calculate pvalue using cdf function
        self.pval_i
            .extend(self.xirs_i.iter().map(|x| 1.0 - chisq.cdf(*x)));
    }

    pub fn finish(&mut self) -> XirsResult {
        // println!("calculate 4. raw");
        self.calculate_raw();
        // println!("calculate 5. norm");
        self.calculate_norm();
        // println!("calculate 6. xirs");
        self.calculate_xirs();
        // println!("calculate 7. xirs");
        self.calculate_pval();
        // println!("calculate 8. writing");

        // println!("cm: {}, {:?}..", self.cm_j.len(), &self.rm_i[1..5]);
        // println!("rm: {}, {:?}..", self.rm_i.len(), &self.rm_i[1..5]);
        // println!("rs: {}, {:?}..", self.rs_i.len(), &self.rs_i[1..5]);
        // println!("xs: {}, {:?}..", self.xs_i.len(), &self.xs_i[1..5]);
        // println!("raw: {}, {:?}..", self.raw_i.len(), &self.raw_i[1..5]);
        // println!("afreq: {}, {:?}..", self.afrq_i.len(), &self.afrq_i[1..5]);
        // println!("xris: {}, {:?}..", self.xirs_i.len(), &self.xirs_i[1..5]);

        let m = self.nsites;
        let mut chr_id = Vec::<u32>::with_capacity(m);
        let mut chr_pos = Vec::<u32>::with_capacity(m);
        let gw_pos = self.pos_i.clone();

        for gwpos in gw_pos.iter() {
            let (chid, _, chrpos) = self.ibd.get_ginfo().to_chr_pos(*gwpos);
            chr_id.push(chid as u32);
            chr_pos.push(chrpos);
        }
        let xirs = self.xirs_i.clone();
        let pval = self.pval_i.clone();
        let raw = self.raw_i.clone();

        XirsResult {
            chr_id,
            chr_pos,
            gw_pos,
            raw,
            xirs,
            pval,
        }
    }
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
            ("Raw", take(&mut self.raw).into_arrow_array()),
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
