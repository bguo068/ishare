use std::ops::Div;

use crate::container::intervals::Intervals;
use crate::container::intervaltree::IntervalTree;
use crate::genome::GenomeInfo;
use crate::share::ibd::coverage::CovCounter;
use crate::stat::xirs::XirsResult;
pub struct IbdPeakFinder<'a> {
    ginfo: &'a GenomeInfo,
    pos: Vec<u32>,
    cov: Vec<usize>,
}

impl<'a> IbdPeakFinder<'a> {
    pub fn new(cov_counter: &'a CovCounter, ginfo: &'a GenomeInfo) -> Self {
        let cov = cov_counter.get_counts().to_owned();
        let pos: Vec<u32> = cov_counter
            .get_intervals()
            .into_iter()
            .map(|r| r.start)
            .collect();
        Self { ginfo, cov, pos }
    }

    pub fn find(&self) -> Intervals<u32> {
        let num_chroms = self.ginfo.chromsize.len();
        let mut intervals_vec = vec![Intervals::new(); num_chroms];
        let mut cov_vec = vec![vec![]; num_chroms];

        // build interval/coverage for each chromosome
        for (gw_pos, cov) in self.pos.iter().zip(self.cov.iter()) {
            let (chrid, _chrname, _chr_pos) = self.ginfo.to_chr_pos(*gw_pos);
            let chr_itvrl = &mut intervals_vec[chrid];
            chr_itvrl.push(*gw_pos..(*gw_pos + 1));
            cov_vec[chrid].push(*cov);
        }

        // expand the ends of intervals so it is easier for merging
        for (chrid, chr_itvls) in intervals_vec.iter_mut().enumerate() {
            let gw_chr_start = self.ginfo.gwstarts[chrid];
            let chr_len = self.ginfo.chromsize[chrid];
            let gw_chr_end = gw_chr_start + chr_len;
            let it = chr_itvls
                .iter()
                .zip(chr_itvls.iter().skip(1))
                .map(|(r1, r2)| (r1.start, r2.start));
            let mut intvs = Intervals::new();
            intvs.extend_from_iter(it);

            let last_start = chr_itvls.iter().last().unwrap().start;
            intvs.push(last_start..gw_chr_end);
            *chr_itvls = intvs;
        }

        // find threshold for each chromosome
        let mut thres_vec = vec![0usize; num_chroms];
        let mut median_vec = vec![0usize; num_chroms];
        for (chrid, chrcov) in cov_vec.iter_mut().enumerate() {
            chrcov.sort();
            // trimming 5% on each sides
            let s = chrcov.len() / 20;
            let e = chrcov.len() - chrcov.len() / 20;

            // get median
            let median = chrcov[(s + e) / 2];
            median_vec[chrid] = median;

            // get mean and std
            let mean = chrcov[s..e].iter().sum::<usize>() as f64 / (e - s) as f64;
            let std = &chrcov[s..e]
                .iter()
                .map(|x| (*x as f64 - mean) * (*x as f64 - mean))
                .sum::<f64>()
                .div((e - s) as f64)
                .sqrt();

            let thres = (mean + 2.0 * std) as usize;
            thres_vec[chrid] = thres;
        }

        // get peak
        let mut peaks = Intervals::new();

        let mut tree = IntervalTree::new(100);
        for chrid in 0..num_chroms {
            let median = median_vec[chrid];
            let thres = thres_vec[chrid];

            let mut core = Intervals::new();
            let mut extension = Intervals::new();
            for (r, cov) in intervals_vec[chrid].iter().zip(cov_vec[chrid].iter()) {
                if *cov > median {
                    extension.push(r.to_owned());
                }
                if *cov > thres {
                    core.push(r.to_owned());
                }
            }

            core.merge();
            extension.merge();

            // any extension regon (merged) that overlap with any peak region (merged)
            // is the extended peak
            tree.clear_and_fill_with_iter(core.iter().map(|r| (r.to_owned(), ())));

            for ext in extension.iter() {
                if tree.query(ext.to_owned()).count() > 0 {
                    peaks.push(ext.to_owned());
                }
            }
        }
        peaks
    }
}

pub fn filter_peaks(peaks: &Intervals<u32>, xirs: &XirsResult) -> Intervals<u32> {
    struct RankedPvalue {
        oid: usize,
        pval: f64,
        rank: usize,
    }

    // adjust p values
    let n = xirs.pval.len();
    let mut pvalues = Vec::<RankedPvalue>::with_capacity(n);
    for (i, pval) in xirs.pval.iter().enumerate() {
        pvalues.push(RankedPvalue {
            oid: i,
            pval: *pval,
            rank: 0,
        })
    }

    pvalues.sort_by(|a, b| a.pval.partial_cmp(&b.pval).unwrap());

    pvalues.iter_mut().enumerate().for_each(|(i, x)| {
        x.rank = i;
        x.pval = x.pval * n as f64 / x.rank as f64;
    });

    pvalues.sort_by_key(|x| x.oid);

    // hits
    let tree = IntervalTree::from_iter(
        pvalues
            .iter()
            .zip(xirs.gw_pos.iter())
            .map(|(rp, gwpos)| (rp.pval < 0.05, gwpos))
            .filter(|x| x.0)
            .map(|x| (*x.1..(*x.1 + 1), ())),
    );

    // peak regions contains hits
    let mut filt_peaks = Intervals::new();

    peaks
        .iter()
        .filter(|r| tree.query((*r).clone()).count() > 0)
        .for_each(|r| filt_peaks.push(r.to_owned()));

    filt_peaks
}
