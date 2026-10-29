use super::utils::{bcf_reader, create_parent_dir};
use super::{BcftkError, Result};
use error_stack::*;

use super::Commands;
use bitvec::prelude::*;
use itertools::{EitherOrBoth, Itertools};
use rust_htslib::bcf::{record::GenotypeAllele, Read};
use std::{io::BufWriter, path::Path};

#[derive(Default)]
pub struct GenotypeMissingnessMatrix {
    data: BitVec<usize, Msb0>,
    nrow: usize,
    ncol: usize,
    sample_mask: BitVec<usize, Msb0>,
    site_mask: BitVec<usize, Msb0>,
    num_sample_used: usize,
    num_site_used: usize,
    sample_nonmiss_count: Vec<(usize, usize)>,
    site_nonmiss_count: Vec<(usize, usize)>,
    // information for recovery site and samples
    contigs: NameVec,
    samples: NameVec,
    site_rids: Vec<u32>,
    site_poss: Vec<i64>,
}

#[derive(Default)]
struct NameVec {
    data: Vec<u8>,
    name_starts: Vec<usize>,
}

impl NameVec {
    fn add(&mut self, name: &[u8]) {
        self.name_starts.push(self.data.len());
        self.data.extend(name.iter());
    }
    fn get(&self, idx: usize) -> Result<&str> {
        let start_byte = self.name_starts[idx];
        let mut end_byte = self.data.len();
        if idx < self.name_starts.len() - 1 {
            end_byte = self.name_starts[idx + 1];
        }
        std::str::from_utf8(&self.data[start_byte..end_byte])
            .change_context(BcftkError::Input)
            .attach("cannot convert name bytes to str")
    }
    fn len(&self) -> usize {
        self.name_starts.len()
    }
}

impl GenotypeMissingnessMatrix {
    pub fn new(input: impl AsRef<Path>) -> Result<Self> {
        let mut m = GenotypeMissingnessMatrix::default();

        let mut reader = bcf_reader(Some(&input.as_ref().to_path_buf()))?;

        let header = reader.header();

        header.samples().iter().for_each(|s| m.samples.add(s));
        for i in 0..header.contig_count() {
            let contig_name = header
                .rid2name(i)
                .change_context(BcftkError::Input)
                .attach("cannot get contig name by id")?;
            m.contigs.add(contig_name);
        }

        let nsam = m.samples.len();
        m.ncol = nsam;
        m.sample_mask = BitVec::<usize, Msb0>::with_capacity(nsam);
        m.sample_mask.resize(nsam, true);

        let mut record = reader.empty_record();

        while let Some(res) = reader.read(&mut record) {
            res.change_context(BcftkError::Input)
                .attach("cannot read a bcf record")?;
            let genotypes = record
                .genotypes()
                .change_context(BcftkError::Input)
                .attach("cannot parse genotypes from a bcf record")?;
            for isample in 0..nsam {
                let missing = genotypes.get(isample).iter().any(|allele| {
                    matches!(
                        allele,
                        GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing
                    )
                });
                m.data.push(!missing);
            }
            m.nrow += 1;
            m.site_rids.push(
                record
                    .rid()
                    .ok_or(BcftkError::Input)
                    .attach("bcf record has empty rid")?,
            );
            m.site_poss.push(record.pos())
        }

        m.site_mask = BitVec::<usize, Msb0>::with_capacity(m.nrow);
        m.site_mask.resize(m.nrow, true);

        m.count_nonmissing();

        Ok(m)
    }

    fn count_nonmissing(&mut self) {
        // borrow the nonmissing count vectors
        let mut sample_nmc = vec![];
        let mut site_nmc = vec![];
        std::mem::swap(&mut sample_nmc, &mut self.sample_nonmiss_count);
        std::mem::swap(&mut site_nmc, &mut self.site_nonmiss_count);
        let mut num_sample_used = 0;
        let mut num_site_used = 0;

        // init
        sample_nmc.clear();
        site_nmc.clear();
        for isample in self.sample_mask.iter_ones() {
            sample_nmc.push((isample, 0));
            num_sample_used += 1;
        }
        for isite in self.site_mask.iter_ones() {
            site_nmc.push((isite, 0));
            num_site_used += 1;
        }
        self.num_sample_used = num_sample_used;
        self.num_site_used = num_site_used;

        // do the counting
        for (isite, isite_nmr) in site_nmc.iter_mut() {
            let site_miss_bit_slice = &self.data[(*isite * self.ncol)..((*isite + 1) * self.ncol)];

            sample_nmc
                .iter_mut()
                .merge_join_by(
                    site_miss_bit_slice.iter().enumerate(),
                    |(isample, _isample_nmr), (ibit, _not_missing)| isample.cmp(ibit),
                )
                .for_each(|res| match res {
                    EitherOrBoth::Left(_) => {
                        unreachable!("None case should be impossible here")
                    }
                    EitherOrBoth::Right(_) => {} // not used sample
                    EitherOrBoth::Both((_isample, isample_nmr), (_ibit, not_missing)) => {
                        if *not_missing {
                            *isample_nmr += 1;
                            *isite_nmr += 1;
                        }
                    }
                });
        }

        // return the nonmissing count vectors
        std::mem::swap(&mut sample_nmc, &mut self.sample_nonmiss_count);
        std::mem::swap(&mut site_nmc, &mut self.site_nonmiss_count);
    }

    pub fn filter_site(&mut self, min_site_nonmissing_rate: f64) {
        // borrow masks
        let mut site_mask = BitVec::<usize, Msb0>::new();
        std::mem::swap(&mut site_mask, &mut self.site_mask);

        // update masks
        let minreq = self.num_sample_used as f64 * min_site_nonmissing_rate;
        for (isite, site_nmc) in self.site_nonmiss_count.iter() {
            if ((*site_nmc as f64) < minreq) || (self.num_sample_used == 0) {
                site_mask.set(*isite, false);
            }
        }

        // return mask
        std::mem::swap(&mut site_mask, &mut self.site_mask);

        // recount nonmissing
        self.count_nonmissing();
    }

    pub fn filter_sample(&mut self, min_sample_nonmissing_rate: f64) {
        // borrow masks
        let mut sample_mask = BitVec::<usize, Msb0>::new();
        std::mem::swap(&mut sample_mask, &mut self.sample_mask);

        // update masks
        let min_req = self.num_site_used as f64 * min_sample_nonmissing_rate;
        for (isample, sample_nmc) in self.sample_nonmiss_count.iter() {
            if ((*sample_nmc as f64) < min_req) || (self.num_site_used == 0) {
                sample_mask.set(*isample, false);
            }
        }

        // return mask
        std::mem::swap(&mut sample_mask, &mut self.sample_mask);

        // recount nonmissing
        self.count_nonmissing();
    }

    pub fn write_used_sites(&self, p: impl AsRef<Path>) -> Result<()> {
        create_parent_dir(&p)?;
        let mut f = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .change_context(BcftkError::Output)
            .attach_with(|| format!("cannot create sites file {}", p.as_ref().display()))?;
        use std::io::Write;

        for ((rid, pos), used) in self
            .site_rids
            .iter()
            .zip(self.site_poss.iter())
            .zip(self.site_mask.iter())
        {
            if *used {
                writeln!(f, "{}:{}", self.contigs.get(*rid as usize)?, *pos + 1)
                    .change_context(BcftkError::Output)?;
            }
        }
        f.flush().change_context(BcftkError::Output)?;
        Ok(())
    }

    pub fn write_used_samples(&self, p: impl AsRef<Path>) -> Result<()> {
        create_parent_dir(&p)?;
        let mut f = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .change_context(BcftkError::Output)
            .attach_with(|| format!("cannot create samples file {}", p.as_ref().display()))?;
        use std::io::Write;
        for isample in self.sample_mask.iter_ones() {
            writeln!(f, "{}", self.samples.get(isample)?).change_context(BcftkError::Output)?;
        }
        f.flush().change_context(BcftkError::Output)?;
        Ok(())
    }
}

#[derive(Debug, PartialEq)]
enum NonMissingRateFilter {
    MinVmissRate(f64),
    MinImissRate(f64),
}

fn parse_filter(filter_strings: &[String]) -> Result<Vec<NonMissingRateFilter>> {
    let mut filters = vec![];
    for filter_str in filter_strings {
        let mut subfields = filter_str.split(":");
        let rate: f64 = subfields
            .next()
            .ok_or(BcftkError::Input)
            .attach_with(|| format!("cannot get rate field from filter: {filter_str}"))?
            .parse()
            .change_context(BcftkError::Input)
            .attach_with(|| format!("cannot parse rate field from filter: {filter_str}"))?;
        let typ = subfields
            .next()
            .ok_or(BcftkError::Input)
            .attach_with(|| format!("cannot get type field from filter: {filter_str}"))?;

        let filter = match typ {
            "v" => NonMissingRateFilter::MinVmissRate(rate),
            "i" => NonMissingRateFilter::MinImissRate(rate),
            _ => bail!(BcftkError::Input.into_report().attach(format!(
                "filter type should be `v` or `i` but got `{typ}` in filter: {filter_str}"
            ))),
        };
        filters.push(filter);
    }
    Ok(filters)
}

pub fn main_filtmiss(args: &Commands) -> Result<()> {
    if let Commands::Filtmiss {
        input,
        filter,
        output_prefix,
    } = args
    {
        let filters = parse_filter(filter)?;
        let mut gm = GenotypeMissingnessMatrix::new(input)?;

        for filter in filters {
            match filter {
                NonMissingRateFilter::MinVmissRate(min_vmiss_rate) => {
                    gm.filter_site(min_vmiss_rate);
                }
                NonMissingRateFilter::MinImissRate(min_imiss_rate) => {
                    gm.filter_sample(min_imiss_rate);
                }
            }
        }

        gm.write_used_sites(format!("{output_prefix}_sites.txt"))?;
        gm.write_used_samples(format!("{output_prefix}_samples.txt"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// build a matrix from row-major (site x sample) nonmissing indicators
    fn from_nonmissing(nonmissing: &[&[bool]]) -> GenotypeMissingnessMatrix {
        let nsite = nonmissing.len();
        let nsam = nonmissing[0].len();
        let mut m = GenotypeMissingnessMatrix {
            ncol: nsam,
            nrow: nsite,
            ..Default::default()
        };
        m.sample_mask.resize(nsam, true);
        m.site_mask.resize(nsite, true);
        for row in nonmissing {
            m.data.extend(row.iter().copied());
        }
        m.contigs.add(b"chr1");
        (0..nsam).for_each(|isample| m.samples.add(format!("s{isample}").as_bytes()));
        m.site_rids.resize(nsite, 0);
        (0..nsite).for_each(|isite| m.site_poss.push(isite as i64 * 10));
        m.count_nonmissing();
        m
    }

    #[test]
    fn test_parse_filter() {
        let fs = parse_filter(&["0.3:v".to_owned(), "0.4:i".to_owned()]).unwrap();
        assert_eq!(
            fs,
            vec![
                NonMissingRateFilter::MinVmissRate(0.3),
                NonMissingRateFilter::MinImissRate(0.4)
            ]
        );
        assert!(parse_filter(&["0.3:x".to_owned()]).is_err());
        assert!(parse_filter(&["0.3".to_owned()]).is_err());
        assert!(parse_filter(&["abc:v".to_owned()]).is_err());
    }

    #[test]
    fn test_filter_site_then_sample() {
        let (t, f) = (true, false);
        // 4 sites x 3 samples
        let mut m = from_nonmissing(&[
            &[t, t, f], // site0: 2/3
            &[t, t, t], // site1: 3/3
            &[f, f, t], // site2: 1/3
            &[t, f, f], // site3: 1/3
        ]);
        assert_eq!(m.site_nonmiss_count, vec![(0, 2), (1, 3), (2, 1), (3, 1)]);
        assert_eq!(m.sample_nonmiss_count, vec![(0, 3), (1, 2), (2, 2)]);

        // keep sites with nonmissing rate >= 0.5: site0, site1
        m.filter_site(0.5);
        assert_eq!(m.site_mask.iter_ones().collect::<Vec<_>>(), vec![0, 1]);
        assert_eq!(m.sample_nonmiss_count, vec![(0, 2), (1, 2), (2, 1)]);

        // keep samples with nonmissing rate >= 0.8 among the kept sites
        m.filter_sample(0.8);
        assert_eq!(m.sample_mask.iter_ones().collect::<Vec<_>>(), vec![0, 1]);
        assert_eq!(m.site_nonmiss_count, vec![(0, 2), (1, 2)]);

        let tmp = tempfile::tempdir().unwrap();
        let sites = tmp.path().join("sub/out_sites.txt");
        let samples = tmp.path().join("sub/out_samples.txt");
        m.write_used_sites(&sites).unwrap();
        m.write_used_samples(&samples).unwrap();
        assert_eq!(std::fs::read_to_string(sites).unwrap(), "chr1:1\nchr1:11\n");
        assert_eq!(std::fs::read_to_string(samples).unwrap(), "s0\ns1\n");
    }
}
