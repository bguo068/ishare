use super::utils::{bcf_reader, create_parent_dir};
use super::{BcftkError, Result};
use error_stack::*;

use super::Commands;
use rust_htslib::bcf::Read;
use std::io::Write;
use std::path::Path;

const BINS: [f64; 12] = [
    -0.0001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.5, 0.5001,
];

pub struct FwsCalculator {
    pop_het_sum: Vec<f64>,      // length: nbins
    pop_het_cnt: Vec<usize>,    // length: nbins
    sample_het_sum: Vec<f64>,   // length: samples x nbins
    sample_het_cnt: Vec<usize>, // length: samples x nbins
    nsamples: usize,
}

impl FwsCalculator {
    pub fn new(nsamples: usize) -> Self {
        Self {
            pop_het_sum: vec![0.0; BINS.len() - 1],
            pop_het_cnt: vec![0; BINS.len() - 1],
            sample_het_sum: vec![0.0; nsamples * (BINS.len() - 1)],
            sample_het_cnt: vec![0; nsamples * (BINS.len() - 1)],
            nsamples,
        }
    }

    pub fn update_stats_for_a_site(
        &mut self,
        site_ad: &[u32],
        nsamples: usize,
        counter: &mut Vec<usize>,
    ) {
        let nalleles = site_ad.len() / nsamples;
        let nbins = BINS.len() - 1;

        // count ad across sample for a given size
        counter.clear();
        counter.resize(nalleles, 0);
        let mut total = 0;
        for (iallele, ad) in site_ad.iter().enumerate() {
            counter[iallele % nalleles] += *ad as usize;
            total += *ad;
        }
        counter.sort_by_key(|c| -(*c as isize)); // descending

        // calculate maf
        let maf = if nalleles < 2 {
            0.0
        } else {
            counter[1] as f64 / total as f64
        };
        let which = BINS.partition_point(|x| *x <= maf) - 1;

        // calculate population hets
        let sq_sum: usize = counter.iter().map(|c| *c * *c).sum();
        let sum_sq = total as f64 * total as f64;
        if sum_sq > 0.0 {
            let pop_het = 1.0 - sq_sum as f64 / sum_sq;
            self.pop_het_sum[which] += pop_het;
            self.pop_het_cnt[which] += 1;
        }

        // get sample hets
        for (isample, chunk) in site_ad.chunks(nalleles).enumerate() {
            let sq_sum = chunk.iter().map(|&c| c as f64 * c as f64).sum::<f64>();
            let sum = chunk.iter().sum::<u32>();
            if sum > 0 {
                let sample_het = 1.0 - sq_sum / sum as f64 / sum as f64;
                self.sample_het_sum[isample * nbins + which] += sample_het;
                self.sample_het_cnt[isample * nbins + which] += 1;
            }
        }
    }

    pub fn calculate_fws(&mut self) -> Vec<f64> {
        let nbins = BINS.len() - 1;
        let mut fws_vec = vec![f64::NAN; self.nsamples];

        for (sum, cnt) in self.pop_het_sum.iter_mut().zip(self.pop_het_cnt.iter_mut()) {
            if *cnt == 0 {
                *sum = f64::NAN;
            } else {
                *sum /= *cnt as f64;
                *cnt = 1;
            }
        }

        for (sum, cnt) in self
            .sample_het_sum
            .iter_mut()
            .zip(self.sample_het_cnt.iter_mut())
        {
            if *cnt == 0 {
                *sum = f64::NAN;
            } else {
                *sum /= *cnt as f64;
                *cnt = 1;
            }
        }

        for (isample, sample_hets) in self.sample_het_sum.chunks(nbins).enumerate() {
            let mut nvalid = 0;
            let mut xy_sum = 0f64;
            let mut xx_sum = 0f64;
            for (sample_het, pop_het) in sample_hets.iter().zip(self.pop_het_sum.iter()) {
                if sample_het.is_nan() || pop_het.is_nan() {
                    continue;
                }
                xx_sum += *pop_het * *pop_het;
                xy_sum += *pop_het * *sample_het;
                nvalid += 1;
            }
            if (xx_sum != 0.0) && (nvalid >= 2) {
                fws_vec[isample] = 1.0 - xy_sum / xx_sum;
            }
        }

        fws_vec
    }
}

pub fn main_calcfws(args: &Commands) -> Result<()> {
    if let Commands::Calcfws { input, output } = args {
        let mut bcfreader = bcf_reader(input.as_ref())?;
        let sample_names_vec: Vec<String> = bcfreader
            .header()
            .samples()
            .iter()
            .map(|s| std::str::from_utf8(s).map(|v| v.to_owned()))
            .collect::<std::result::Result<Vec<_>, _>>()
            .change_context(BcftkError::Input)
            .attach("cannot read sample names from bcf")?;
        let nsam = sample_names_vec.len();
        let mut record = bcfreader.empty_record();
        let mut fws_calculator = FwsCalculator::new(nsam);

        let mut site_ad = vec![0u32; 0];
        let mut counter = vec![];
        let mut nsite_skipped = 0;

        while let Some(res) = bcfreader.read(&mut record) {
            res.change_context(BcftkError::Input)
                .attach("cannot parse bcf record")?;
            site_ad.clear();

            for sample_ad in record
                .format(b"AD")
                .integer()
                .change_context(BcftkError::Input)
                .attach("cannot read FMT/AD field")?
                .iter()
            {
                sample_ad.iter().for_each(|ad| site_ad.push(*ad as u32));
            }
            // skip multiallelic sites
            if site_ad.len() != 2 * nsam {
                nsite_skipped += 1;
                continue;
            }
            fws_calculator.update_stats_for_a_site(&site_ad, nsam, &mut counter);
        }
        if nsite_skipped > 0 {
            eprintln!("WARN: {nsite_skipped} sites are skipped as they are not biallelic");
        }
        let fws_vec = fws_calculator.calculate_fws();
        let outpath = output
            .as_ref()
            .map(|p| p.as_path())
            .unwrap_or(Path::new("out_fws.tsv"));
        create_parent_dir(outpath)?;
        let mut outfile = std::fs::File::create(outpath)
            .map(std::io::BufWriter::new)
            .change_context(BcftkError::Output)
            .attach_with(|| format!("cannot create output file {}", outpath.display()))?;

        for (sample_name, fws) in sample_names_vec.iter().zip(fws_vec.iter()) {
            writeln!(outfile, "{sample_name}\t{fws}").change_context(BcftkError::Output)?;
        }
        outfile.flush().change_context(BcftkError::Output)?;
        eprintln!("INFO: wrote to file:  {}", outpath.display());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::PathBuf;

    fn read_fws(fws_path: impl AsRef<Path>) -> HashMap<String, f64> {
        let mut res = HashMap::new();
        for line in std::fs::read_to_string(fws_path).unwrap().lines() {
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let sample_name = fields.next().unwrap().to_owned();
            let fws_val: f64 = fields.next().unwrap().parse().unwrap();
            res.insert(sample_name, fws_val);
        }
        res
    }

    #[test]
    fn test_fws_clonal_vs_mixed() {
        // sample 0 is clonal (all reads support one allele per site);
        // sample 1 is mixed (reads split evenly at each site)
        let nsam = 2;
        let mut calc = FwsCalculator::new(nsam);
        let mut counter = vec![];
        for (a, b) in [(10u32, 0u32), (0, 10), (10, 0), (0, 10)] {
            // population maf varies a little across sites to fill > 1 bins
            for extra in [0u32, 4] {
                let site_ad = [a, b, 5 + extra, 5];
                calc.update_stats_for_a_site(&site_ad, nsam, &mut counter);
            }
        }
        let fws = calc.calculate_fws();
        assert_eq!(fws.len(), 2);
        assert!(fws[0] > 0.9, "clonal sample fws={}", fws[0]);
        assert!(fws[1] < fws[0], "mixed sample fws={}", fws[1]);
    }

    #[test]
    fn test_calcfws_matches_moimix() {
        let testdata = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../testdata/dir003");
        let tmp = tempfile::tempdir().unwrap();
        let output = tmp.path().join("sub/fws_calc.tsv");

        main_calcfws(&Commands::Calcfws {
            input: Some(testdata.join("prefilt.vcf.gz")),
            output: Some(output.clone()),
        })
        .unwrap();

        let fwscal_res = read_fws(&output);
        let moimix_res = read_fws(testdata.join("moimix_fws.tsv"));
        assert_eq!(fwscal_res.len(), moimix_res.len());

        for (sample, fwscal_val) in fwscal_res {
            let moimix_val = moimix_res[sample.as_str()];
            let rel_diff = (moimix_val - fwscal_val).abs() / (moimix_val + fwscal_val);
            assert!(
                rel_diff < 0.001,
                "inconsistent result for {sample}. fwscal: {fwscal_val}, moimix: {moimix_val}"
            );
        }
    }
}
