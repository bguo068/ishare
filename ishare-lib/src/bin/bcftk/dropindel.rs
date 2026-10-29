use super::utils::{bcf_reader, create_parent_dir, parse_bcf_output_format};
use super::{BcftkError, Result};
use error_stack::*;

use super::Commands;
use rust_htslib::bcf::{self, record::GenotypeAllele, Read};

/// Helper function determine if an alternative allele is a snp relative to the
/// ref allele.
///
/// The first argument are slice of alleles (from `record.alleles`) which
/// include both ref and all alt alleles from this record. the first alt allele
/// corresponding to which = 1. if which = 0 or which >= the number of all
/// alleles, return None; otherwise return Some (bool) with the boolean
/// indicating whether the alt allele is snp relative to ref allele.
///
/// Note: the function assumes allele characters are single-byte.
fn is_alt_allele_snp(alleles: &[&[u8]], which: usize) -> Option<bool> {
    if (which == 0) || (which >= alleles.len()) {
        None
    } else {
        let ref_allele = alleles[0];
        let alt_allele = alleles[which];
        if ref_allele.len() != alt_allele.len() {
            Some(false)
        } else {
            let num_diff_nucleotides = ref_allele
                .iter()
                .zip(alt_allele.iter())
                .filter(|(ref_byte, alt_byte)| **ref_byte != **alt_byte)
                .count();
            Some(num_diff_nucleotides == 1)
        }
    }
}

fn recode_alleles<'a>(alleles: &[&'a [u8]], new_alleles: &mut Vec<&'a [u8]>, keep_allele: &[bool]) {
    new_alleles.clear();
    for (allele, _keep) in alleles
        .iter()
        .zip(keep_allele.iter())
        .filter(|(_, keep)| **keep)
    {
        new_alleles.push(*allele);
    }
}

/// assume the length of keep_allele is the same as that of `record.alleles()`
fn recode_fmt_gt(genotype: &mut [GenotypeAllele], keep_allele: &[bool]) {
    const SENTINEL_GT: i32 = 99999;
    let index_map = {
        let mut index_map = vec![];
        let mut new_idx = 0;
        for keep in keep_allele {
            if *keep {
                index_map.push(new_idx);
                new_idx += 1;
            } else {
                index_map.push(SENTINEL_GT);
            }
        }
        index_map
    };
    // remove gt for alleles not kept
    genotype.iter_mut().for_each(|gta| {
        *gta = match *gta {
            GenotypeAllele::Unphased(orig) if index_map[orig as usize] == SENTINEL_GT => {
                GenotypeAllele::UnphasedMissing
            }
            GenotypeAllele::Unphased(orig) => GenotypeAllele::Unphased(index_map[orig as usize]),
            GenotypeAllele::Phased(orig) if index_map[orig as usize] == SENTINEL_GT => {
                GenotypeAllele::PhasedMissing
            }
            GenotypeAllele::Phased(orig) => GenotypeAllele::Phased(index_map[orig as usize]),
            x => x,
        }
    });
}

fn recode_fmt_ad(ad: &mut Vec<i32>, keep_allele: &[bool]) {
    let nalleles = keep_allele.len();
    const SENTINEL_AD: i32 = i32::MAX;
    for chunk in ad.chunks_mut(nalleles) {
        chunk
            .iter_mut()
            .zip(keep_allele.iter())
            .for_each(|(allele_ad, keep)| {
                if !*keep {
                    *allele_ad = SENTINEL_AD;
                }
            });
    }
    ad.retain(|x| *x != SENTINEL_AD);
}

fn calc_combined_nonsnp_afreq(
    gt: &[GenotypeAllele],
    snp_allele: &[bool],
    allele_counter: &mut Vec<usize>,
) -> f64 {
    allele_counter.clear();
    allele_counter.resize(snp_allele.len(), 0);
    let mut total_count = 0;
    let mut indel_count = 0;
    for allele_idx in gt.iter().filter_map(|s| s.index()) {
        total_count += 1;
        allele_counter[allele_idx as usize] += 1;
    }
    for (c, is_snp) in allele_counter.iter().zip(snp_allele.iter()) {
        if !*is_snp {
            indel_count += *c;
        }
    }
    if total_count == 0 {
        f64::NAN
    } else {
        indel_count as f64 / total_count as f64
    }
}

pub fn main_dropindel(args: &Commands) -> Result<()> {
    if let Commands::Dropindel {
        input,
        output,
        output_format,
        max_total_nonsnp_afreq,
    } = args
    {
        let mut reader = bcf_reader(input.as_ref())?;
        let header_w = bcf::Header::from_template(reader.header());
        let (out_uncompressed, out_fmt) = parse_bcf_output_format(output_format)?;
        let mut writer = match output.as_ref() {
            Some(output) => {
                create_parent_dir(output)?;
                bcf::Writer::from_path(output, &header_w, out_uncompressed, out_fmt)
                    .change_context(BcftkError::Output)
                    .attach_with(|| format!("cannot create bcf writer from {}", output.display()))?
            }
            None => bcf::Writer::from_stdout(&header_w, out_uncompressed, out_fmt)
                .change_context(BcftkError::Output)
                .attach("cannot create bcf writer from stdout")?,
        };
        let mut record = reader.empty_record();
        let mut record_w = writer.empty_record();

        let mut keep_allele = vec![];
        let mut gt: Vec<GenotypeAllele> = vec![];
        let mut ad: Vec<i32> = vec![];
        let mut allele_counter = vec![];
        let mut ploidy: Option<usize> = None;

        while let Some(res) = reader.read(&mut record) {
            res.change_context(BcftkError::Input)
                .attach("cannot read a record from bcf reader")?;
            let mut new_alleles = vec![];
            gt.clear();

            // calculate keep allele vector
            keep_allele.clear();
            keep_allele.push(true);
            let alleles = record.alleles();
            for which_allele in 1..alleles.len() {
                keep_allele.push(
                    is_alt_allele_snp(&alleles, which_allele)
                        .ok_or(BcftkError::Input)
                        .attach("invalid alternative allele index")?,
                );
            }

            // get gt vector
            for sample_gt in record
                .format(b"GT")
                .integer()
                .change_context(BcftkError::Input)
                .attach("cannot interpret FMT/GT as integers")?
                .iter()
            {
                // check ploidy is consistent across samples and sites
                match ploidy {
                    Some(p) if p != sample_gt.len() => {
                        bail!(BcftkError::Input
                            .into_report()
                            .attach(format!("inconsistent ploidy: {p} vs {}", sample_gt.len())))
                    }
                    Some(_) => {}
                    None => ploidy = Some(sample_gt.len()),
                }
                gt.extend(sample_gt.iter().map(|s| GenotypeAllele::from(*s)));
            }

            // check total allele frequency of non snp alleles
            //
            // if total allele frequency is nan or too high skip of this site (filtered out)
            if let Some(max_total_nonsnp_afreq) = max_total_nonsnp_afreq.as_ref() {
                let total_nonsnp_afreq =
                    calc_combined_nonsnp_afreq(&gt, &keep_allele, &mut allele_counter);
                if total_nonsnp_afreq.is_nan() || (total_nonsnp_afreq > *max_total_nonsnp_afreq) {
                    continue;
                }
            }

            // get ad vector
            ad.clear();
            ad.extend(
                record
                    .format(b"AD")
                    .integer()
                    .change_context(BcftkError::Input)
                    .attach("cannot interpret FMT/AD as integers")?
                    .iter()
                    .flat_map(|s| s.iter()),
            );

            // recode
            recode_alleles(&alleles, &mut new_alleles, &keep_allele);
            recode_fmt_ad(&mut ad, &keep_allele);
            recode_fmt_gt(&mut gt, &keep_allele);

            // prepare new record to write
            record_w.set_rid(record.rid());
            record_w.set_pos(record.pos());
            record_w
                .set_alleles(&new_alleles)
                .change_context(BcftkError::Library)
                .attach("cannot set alleles")?;
            record_w
                .push_genotypes(&gt)
                .change_context(BcftkError::Library)
                .attach("cannot push genotypes")?;
            record_w
                .push_format_integer(b"AD", &ad)
                .change_context(BcftkError::Library)
                .attach("cannot push FMT/AD")?;
            for filt in record.filters() {
                record_w
                    .push_filter(&filt)
                    .change_context(BcftkError::Library)
                    .attach("cannot push a filter")?;
            }
            writer
                .write(&record_w)
                .change_context(BcftkError::Output)
                .attach("cannot write a bcf record")?;
            record_w.clear();
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};

    fn create_test_bcf(vcf_path: &Path) {
        // prepare header lines
        let mut header = bcf::Header::new();
        for header_record in [
            r#"##contig=<ID=Pf3D7_01_v3,length=640851>"#,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
            r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
        ] {
            header.push_record(header_record.as_bytes());
        }

        for sample in ["sample1", "sample2"] {
            header.push_sample(sample.as_bytes());
        }

        let mut writer = bcf::Writer::from_path(vcf_path, &header, true, bcf::Format::Vcf).unwrap();
        let mut record = writer.empty_record();

        // record 1: ref A, alt AC (indel) and G (snp)
        record.set_rid(Some(0));
        record.set_pos(100);
        record.set_alleles(&[b"A", b"AC", b"G"]).unwrap();
        let gt = vec![
            GenotypeAllele::Unphased(0),
            GenotypeAllele::Unphased(1),
            GenotypeAllele::Unphased(0),
            GenotypeAllele::Unphased(2),
        ];
        record.push_genotypes(&gt).unwrap();
        let ad = vec![100, 30, 0, 80, 0, 30];
        record.push_format_integer(b"AD", &ad).unwrap();
        writer.write(&record).unwrap();

        // record 2: biallelic snp
        record.set_rid(Some(0));
        record.set_pos(300);
        record.set_alleles(&[b"G", b"C"]).unwrap();
        let gt = vec![
            GenotypeAllele::Unphased(0),
            GenotypeAllele::Unphased(1),
            GenotypeAllele::Unphased(0),
            GenotypeAllele::Unphased(1),
        ];
        record.push_genotypes(&gt).unwrap();
        let ad = vec![100, 30, 80, 20];
        record.push_format_integer(b"AD", &ad).unwrap();
        writer.write(&record).unwrap();
    }

    /// returns number of records and whether all records are snps
    fn check_output_vcf(p: &Path) -> (usize, bool) {
        let mut reader = bcf::Reader::from_path(p).unwrap();
        let mut record = reader.empty_record();
        let mut nrec = 0;
        let mut all_snp = true;
        while let Some(res) = reader.read(&mut record) {
            res.unwrap();
            nrec += 1;
            all_snp &= unsafe { rust_htslib::htslib::bcf_is_snp(record.inner) == 1 };
        }
        (nrec, all_snp)
    }

    fn run(input: &Path, output: &Path, max_total_nonsnp_afreq: Option<f64>) {
        main_dropindel(&Commands::Dropindel {
            input: Some(PathBuf::from(input)),
            output: Some(PathBuf::from(output)),
            output_format: "v".to_owned(),
            max_total_nonsnp_afreq,
        })
        .unwrap();
    }

    #[test]
    fn test_is_alt_allele_snp() {
        let alleles: [&[u8]; 4] = [b"A", b"AC", b"G", b"T"];
        assert_eq!(is_alt_allele_snp(&alleles, 0), None);
        assert_eq!(is_alt_allele_snp(&alleles, 1), Some(false));
        assert_eq!(is_alt_allele_snp(&alleles, 2), Some(true));
        assert_eq!(is_alt_allele_snp(&alleles, 4), None);
        let alleles: [&[u8]; 3] = [b"AC", b"GT", b"AT"];
        assert_eq!(is_alt_allele_snp(&alleles, 1), Some(false)); // mnp
        assert_eq!(is_alt_allele_snp(&alleles, 2), Some(true));
    }

    #[test]
    fn test_recode_gt_and_ad() {
        let keep = [true, false, true];
        let mut gt = vec![
            GenotypeAllele::Unphased(0),
            GenotypeAllele::Unphased(1),
            GenotypeAllele::Phased(2),
            GenotypeAllele::UnphasedMissing,
        ];
        recode_fmt_gt(&mut gt, &keep);
        assert_eq!(
            gt,
            vec![
                GenotypeAllele::Unphased(0),
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::Phased(1),
                GenotypeAllele::UnphasedMissing,
            ]
        );
        let mut ad = vec![100, 30, 0, 80, 0, 30];
        recode_fmt_ad(&mut ad, &keep);
        assert_eq!(ad, vec![100, 0, 80, 30]);
    }

    #[test]
    fn test_dropindel() {
        let dir = tempfile::tempdir().unwrap();
        let vcf_input_path = dir.path().join("tmp_input.vcf");
        let vcf_output_path = dir.path().join("sub/tmp_output.vcf");
        create_test_bcf(&vcf_input_path);

        // indel allele dropped but record kept
        run(&vcf_input_path, &vcf_output_path, None);
        let (nrec, all_snp) = check_output_vcf(&vcf_output_path);
        assert_eq!(nrec, 2);
        assert!(all_snp);

        // record 1 has indel allele freq 0.25 > 0.1, so it is removed
        run(&vcf_input_path, &vcf_output_path, Some(0.1));
        let (nrec, all_snp) = check_output_vcf(&vcf_output_path);
        assert_eq!(nrec, 1);
        assert!(all_snp);
    }
}
