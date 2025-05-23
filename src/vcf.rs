use ahash::AHashSet;

use crate::{
    container::intervals::Intervals,
    genome::GenomeInfo,
    genotype::{common::*, rare::*},
    gmap::GeneticMap,
    indiv::Individuals,
    site::*,
};
use rust_htslib::bcf::{record::GenotypeAllele, IndexedReader, Read, Reader};
use std::path::Path;

pub fn read_vcf(
    target_samples: &AHashSet<String>,
    gconfig: &GenomeInfo,
    vcf_path: impl AsRef<Path>,
    max_maf: f64,
    region: Option<(u32, u64, Option<u64>)>,
) -> (Sites, Individuals, GenotypeRecords) {
    // let  vcf_fn =  "/local/chib/TOPMedFull/data/org/87668/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.10b/phased/freeze.10b.chr22.pass_only.phased.bcf";
    // let vcf_fn = "./testdata/test.bcf";
    let mut bcf = IndexedReader::from_path(vcf_path).unwrap();

    if let Some((rid, start, end_opt)) = region {
        let chrname = &gconfig.chromnames[rid as usize];
        let rid2 = bcf.header().name2rid(chrname.as_bytes()).unwrap();
        bcf.fetch(rid2, start, end_opt).unwrap();
    }
    bcf.set_threads(3).unwrap();

    let header = bcf.header().clone();
    let nsam = header.sample_count();

    // Output:
    let mut records = Vec::new();
    let mut sites: Sites = Sites::new();
    // Individuals
    let it = header
        .samples()
        .into_iter()
        .map(|x| std::str::from_utf8(x).unwrap());

    // calculate sample_mask
    let mut sample_mask = vec![true; header.sample_count() as usize];
    if target_samples.len() > 0 {
        for (i, s) in it.clone().enumerate() {
            if !target_samples.contains(s) {
                sample_mask[i] = false;
            }
        }
    }

    // only sameples in targets
    let individuals = Individuals::from_str_iter(
        it.zip(sample_mask.iter())
            .filter(|(_s, yes)| **yes)
            .map(|(s, _yes)| s),
    );
    let sel_nsam = individuals.v().len();

    let ac_thres = (sel_nsam as f64 * max_maf * 2.0).floor() as usize;

    let mut ab = AlleleBuffer::new();
    let mut rid_last = None;
    let mut pos_last: Option<u32> = None;
    let mut gt = Vec::<u8>::new();

    let mut allele_counts = Vec::new();
    let mut allele_is_rare = Vec::<bool>::new();
    for record_result in bcf.records() {
        let record = record_result.unwrap();

        let alleles = record.alleles();

        // TODO: use more strict to define variant type
        // check allele type and skip non-SNP
        let mut is_snp_type = true;
        for a in alleles.iter() {
            let a = *a;
            if (a.len() > 1) || (a == b"*") {
                is_snp_type = false;
            }
        }
        if !is_snp_type {
            continue;
        }

        // chromosome id
        let chrid = {
            let rid = record.rid().unwrap();
            let chrname = header.rid2name(rid).unwrap();
            let chrname = std::str::from_utf8(chrname).unwrap();
            gconfig.idx[chrname]
        };

        // genome-wide pos
        let pos = { record.pos() as u32 + gconfig.gwstarts[chrid] as u32 };

        let mut is_new_pos = true;
        if pos_last.is_some() {
            let rid_last = rid_last.unwrap();
            let pos_last = pos_last.unwrap();
            if rid_last == chrid {
                assert!(pos_last <= pos);
            }
            is_new_pos = !((pos_last == pos) && (rid_last == chrid));
        }

        let start_allele: u8;
        if is_new_pos {
            // output information from old pos
            if pos_last.is_some() {
                let buffers = (
                    &mut sites,
                    &mut records,
                    &mut ab,
                    &mut allele_counts,
                    &mut allele_is_rare,
                );
                let pos_last = pos_last.unwrap();
                output_rare_allele_records(buffers, &gt, pos_last, ac_thres);
            }

            ab.clear();
            start_allele = 0;
            ab.push(alleles[0]);
            for allele in alleles.iter().skip(1) {
                ab.push(b" ");
                ab.push_to_data_only(allele);
            }
        } else {
            // Minus 1 as only non ref are added
            start_allele = ab.len() as u8 - 1;
            // assert different lines of the same position have a consistent REF allele
            let alleles = record.alleles();
            assert_eq!(
                ab.get(0),
                alleles[0],
                "REF alleles are consistent for variants with the same position"
            );
            // skip the ref allele as it has been added
            for allele in alleles.iter().skip(1) {
                ab.push(b" ");
                ab.push_to_data_only(allele);
            }
        }
        let fmt = record.format(b"GT");
        let bcf_fmt = fmt.inner();
        // about bcf genotype encoding
        // see https://github.com/samtools/htslib/blob/99415e2a2ce26bdbf4e910954330ea769de2c3f0/htslib/vcf.h#L156
        // https://samtools.github.io/hts-specs/VCFv4.2.pdf
        // assert_eq!(bcf_fmt.type_, 1);
        assert_eq!(bcf_fmt.n, 2);
        assert_eq!(bcf_fmt.p_len, 2 * nsam);

        let raw_gt_bytes = unsafe { std::slice::from_raw_parts(bcf_fmt.p, bcf_fmt.p_len as usize) };

        // make a vector allele idx
        let mut hap_mask = vec![];
        for sm in sample_mask.iter() {
            hap_mask.push(*sm);
            hap_mask.push(*sm);
        }

        let gt_iter = raw_gt_bytes
            .chunks(bcf_fmt.type_ as usize) // each chunk is an integer
            .zip(hap_mask.iter())
            // skip 'no' haplotypes
            .filter(|(_, &m)| m)
            .map(|(x, _)| {
                // bytes to integr
                let mut y: u32 = 0;
                for (ibyte, byte) in x.iter().enumerate() {
                    y |= (*byte as u32) << (ibyte * 8);
                }
                // interger to genotype
                let mut i = match (y as i32).into() {
                    GenotypeAllele::Unphased(i) => i as u8,
                    GenotypeAllele::Phased(i) => i as u8,
                    _ => {
                        panic!("Currently GT should not be missing!")
                    }
                };
                if i != 0 {
                    // rebase
                    i += start_allele;
                }
                i
            });

        if is_new_pos {
            // read genotype
            gt.clear();
            gt.extend(gt_iter);
        } else {
            // update genotype
            for (a, b) in gt.iter_mut().zip(gt_iter) {
                if b != 0 {
                    *a = b;
                }
            }
        }

        // save to ab_last in case some multiallelic sites are normalized into multiple lines of biallelic variant.
        // ab_last allows to refer alleles the previous lines that has the same position
        rid_last = Some(chrid);
        pos_last = Some(pos);
    }

    // output the last record:
    if pos_last.is_some() {
        let pos_last = pos_last.unwrap();
        let buffers = (
            &mut sites,
            &mut records,
            &mut ab,
            &mut allele_counts,
            &mut allele_is_rare,
        );
        output_rare_allele_records(buffers, &gt, pos_last, ac_thres);
    }

    let gtrec = GenotypeRecords::new(records, 1);

    (sites, individuals, gtrec)
}

pub fn read_vcf_for_genotype_matrix(
    target_samples: &AHashSet<String>,
    gconfig: &GenomeInfo,
    vcf_path: impl AsRef<Path>,
    min_maf: f64,
    region: Option<(u32, u64, Option<u64>)>,
) -> (Sites, Individuals, GenotypeMatrix) {
    use rust_htslib::bcf::{record::GenotypeAllele, IndexedReader, Read};
    // let  vcf_fn =  "/local/chib/TOPMedFull/data/org/87668/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.10b/phased/freeze.10b.chr22.pass_only.phased.bcf";
    // let vcf_fn = "./testdata/test.bcf";
    let mut bcf = IndexedReader::from_path(vcf_path).unwrap();

    if let Some((rid, start, end_opt)) = region {
        let chrname = &gconfig.chromnames[rid as usize];
        let rid2 = bcf.header().name2rid(chrname.as_bytes()).unwrap();
        bcf.fetch(rid2, start, end_opt).unwrap();
    }
    bcf.set_threads(3).unwrap();

    let header = bcf.header().clone();
    let nsam = header.sample_count();

    // Output:
    let mut sites: Sites = Sites::new();
    // Individuals
    let it = header
        .samples()
        .into_iter()
        .map(|x| std::str::from_utf8(x).unwrap());

    // calculate sample_mask
    let mut sample_mask = vec![true; header.sample_count() as usize];
    if target_samples.len() > 0 {
        for (i, s) in it.clone().enumerate() {
            if !target_samples.contains(s) {
                sample_mask[i] = false;
            }
        }
    }

    // only sameples in targets
    let individuals = Individuals::from_str_iter(
        it.zip(sample_mask.iter())
            .filter(|(_s, yes)| **yes)
            .map(|(s, _yes)| s),
    );

    let sel_nsam = individuals.v().len();
    let mut gm = GenotypeMatrix::new((sel_nsam as usize) * 2);

    let ac_thres = (sel_nsam as f64 * min_maf * 2.0).floor() as usize;
    let mut rid_last = None;
    let mut pos_last = None;

    // let mut allele_counts = Vec::new();
    for record_result in bcf.records() {
        let record = record_result.unwrap();

        let alleles = record.alleles();

        // TODO: use more strict to define variant type
        // check allele type and skip non-SNP
        let mut is_biallelic_snp_type = true;
        for (i, a) in alleles.iter().enumerate() {
            let a = *a;
            if i > 1 {
                eprintln!("Encoding vcf into genotype matrix only supports biallelic alleles per line for better bit packing;
                    Considering using `bcftools norm -m- ....` to split multiallelic sites into biallelic records. ");
                is_biallelic_snp_type = false;
            }
            if (a.len() > 1) || (a == b"*") {
                is_biallelic_snp_type = false;
                break;
            }
        }
        if !is_biallelic_snp_type {
            continue;
        }

        // chromosome id
        let chrid = {
            let rid = record.rid().unwrap();
            let chrname = header.rid2name(rid).unwrap();
            let chrname = std::str::from_utf8(chrname).unwrap();
            gconfig.idx[chrname]
        };

        // genome-wide pos
        let pos = { record.pos() as u32 + gconfig.gwstarts[chrid] as u32 };

        // check positions are in order
        let mut _is_new_pos = true;
        if pos_last.is_some() {
            let rid_last = rid_last.unwrap();
            let pos_last = pos_last.unwrap();
            if rid_last == chrid {
                assert!(pos_last <= pos, "VCF is not sorted by position");
            }
            _is_new_pos = !((pos_last == pos) && (rid_last == chrid));
        }

        let fmt = record.format(b"GT");
        let bcf_fmt = fmt.inner();
        // assert_eq!(bcf_fmt.type_, 1); // uint8_t
        assert!(bcf_fmt.n == 2);
        assert_eq!(bcf_fmt.p_len, 2 * nsam);

        let raw_gt_slice = unsafe { std::slice::from_raw_parts(bcf_fmt.p, bcf_fmt.p_len as usize) };

        // make a vector allele idx
        let mut hap_mask = vec![];
        for sm in sample_mask.iter() {
            hap_mask.push(*sm);
            hap_mask.push(*sm);
        }
        let gt_iter = raw_gt_slice
            .chunks(bcf_fmt.type_ as usize)
            .zip(hap_mask.iter())
            // skip 'no' haplotypes
            .filter(|(_, &m)| m)
            .map(|(x, _)| {
                // bytes to integr
                let mut y: u32 = 0;
                for (ibyte, byte) in x.iter().enumerate() {
                    y |= (*byte as u32) << (ibyte * 8);
                }
                let i = match (y as i32).into() {
                    GenotypeAllele::Unphased(i) => i as u8,
                    GenotypeAllele::Phased(i) => i as u8,
                    _ => {
                        panic!("Currently GT should not be missing!")
                    }
                };
                i == 1
            });
        let ac: usize = gt_iter.clone().map(|x| if x { 1 } else { 0 }).sum();

        if (ac < ac_thres) || (((nsam * 2) as usize - ac) < ac_thres) {
            // skip rare variants
            continue;
        }
        gm.extend_gt_calls(gt_iter);
        sites.add_site_with_bytes(pos, alleles[0]);
        sites.append_bytes_to_last_allele(b" ");
        sites.append_bytes_to_last_allele(alleles[1]);

        rid_last = Some(chrid);
        pos_last = Some(pos);
    }

    (sites, individuals, gm)
}

type ParamBufferInfo<'a> = (
    &'a mut Sites,
    &'a mut Vec<GenotypeRecord>,
    &'a mut AlleleBuffer,
    &'a mut Vec<usize>,
    &'a mut Vec<bool>,
    // sites: &mut Sites,
    // records: &mut Vec<GenotypeRecord>,
    // ab: &mut AlleleBuffer,
    // allele_counts: &mut Vec<usize>,
    // allele_is_rare: &mut Vec<bool>,
);

fn output_rare_allele_records(buffers: ParamBufferInfo, gt: &[u8], pos: u32, ac_thres: usize) {
    let (sites, records, ab, allele_counts, allele_is_rare) = buffers;
    // Get allele counts
    // init
    allele_counts.clear();
    allele_is_rare.clear();
    for _ in 0..ab.len() {
        allele_counts.push(0);
    }
    // count
    for a in gt.iter() {
        if *a as usize >= ab.len() {
            println!("{:?}", ab);
        }
        allele_counts[*a as usize] += 1;
    }

    // encode the REF and all ALT alleles with non-zero allele counts
    let mut new_encode = 1u8;
    for (i, ac) in allele_counts.iter().enumerate() {
        if *ac < ac_thres {
            allele_is_rare.push(true);
        } else {
            allele_is_rare.push(false);
        }
        assert!(
            i < u8::MAX as usize - 1,
            "allele index should be less than 255"
        );
        // for ALT allele with allele count > 0
        if *ac > 0 {
            ab.set_enc(i as u8, new_encode);
            new_encode += 1;
        }
        // for ALT allele with allele count == 0
        else {
            ab.set_enc(i as u8, u8::MAX); // mark for cleaning up
        }
    }

    if new_encode == 1 {
        // no rare variant in this site, skip adding the variants and sites
        return;
    }

    // Remap allelic encodings and generate rare variants records:
    for (i, old_enc) in gt.iter().enumerate() {
        // skip encoding common variants for genotype records
        if !allele_is_rare[*old_enc as usize] {
            continue;
        }
        // encode rare variants
        let new_enc = ab.get_enc(*old_enc);
        let mut rec = GenotypeRecord::new(0);
        rec.set(pos, i as u32, new_enc);
        records.push(rec);
    }

    // add sites info
    sites.add_site(pos, ab);
}

pub fn read_pos_from_text_file(file_path: &str, ginfo: &GenomeInfo) -> Vec<u32> {
    let mut out = vec![];
    for line in std::fs::read_to_string(file_path)
        .unwrap()
        .trim()
        .split("\n")
    {
        let mut fields = line.split("\t");
        let chrname = fields.next().unwrap();
        let chr_id = ginfo.idx[chrname];
        let chr_pos = fields.next().unwrap().parse().unwrap();
        let gw_pos = ginfo.to_gw_pos(chr_id, chr_pos);
        out.push(gw_pos);
    }
    out
}

pub fn read_pos_from_vcf_file(vcf_path: &str, ginfo: &GenomeInfo) -> Vec<u32> {
    let mut out = vec![];
    let mut bcf = Reader::from_path(vcf_path).unwrap();
    let mut last_rid = u32::MAX;
    let mut chrid = 0usize;
    for record_res in bcf.records() {
        let record = record_res.unwrap();
        let this_rid = record.rid().unwrap();
        if this_rid != last_rid {
            let rname_bstr = record.header().rid2name(this_rid).unwrap();
            let rname_str = std::str::from_utf8(rname_bstr).unwrap();
            chrid = ginfo.idx[rname_str];
            last_rid = this_rid;
        }
        let gw_pos = ginfo.to_gw_pos(chrid, record.pos() as u32);
        out.push(gw_pos);
    }
    out
}

pub fn find_intervals_with_low_snps_density(
    window_size_cm: f32,
    min_snp_per_cm: f32,
    gw_pos: &[u32],
    ginfo: &GenomeInfo,
    gmap: &GeneticMap,
) -> Intervals<u32> {
    let mut low_snp_dens_intervals = Intervals::<u32>::new();
    // find chromosome boundaries
    let mut cm_boundaries = vec![];
    for start in ginfo.gwstarts.iter() {
        cm_boundaries.push(gmap.get_cm(*start));
    }
    cm_boundaries.push(gmap.get_size_cm());

    // for each chromosome, add windows with size of window_size_cm centimorgans
    let mut bp_win_boundaries = vec![];
    let mut win_size_vec_cm = vec![];
    for cms in cm_boundaries.windows(2) {
        let mut s = cms[0];
        let e_max = cms[1];
        while s < e_max {
            let mut e = s + window_size_cm;
            if e > e_max {
                e = e_max;
            }
            bp_win_boundaries.push(gmap.get_bp(s));
            win_size_vec_cm.push(e - s);
            s = e;
        }
    }
    bp_win_boundaries.push(ginfo.get_total_len_bp());

    // window number is 1 smaller than window boundires
    let mut counter = vec![0; bp_win_boundaries.len() - 1];
    for pos in gw_pos.iter() {
        let idx = bp_win_boundaries.partition_point(|x| x <= pos) - 1;
        counter[idx] += 1;
    }
    assert_eq!(counter.len(), win_size_vec_cm.len());
    for ((c, w), win) in counter
        .iter()
        .zip(win_size_vec_cm.iter())
        .zip(bp_win_boundaries.windows(2))
    {
        if *c as f32 / *w < min_snp_per_cm {
            low_snp_dens_intervals.push(win[0]..win[1]);
        }
    }
    low_snp_dens_intervals.merge();
    low_snp_dens_intervals
}

#[cfg(test)]
mod test {
    use crate::{container::intervals::Intervals, genome::Genome};

    use super::find_intervals_with_low_snps_density;

    fn prepare_data() -> (Genome, Vec<u32>, Intervals<u32>) {
        let genome = Genome::new_from_constant_recombination_rate(
            "test",
            &[10_000_000, 10_000_000, 10_000_000],
            &["chr1".into(), "chr2".into(), "chr3".into()],
            1e-8,
        );
        let pos: Vec<u32> = (50..30_000_000)
            .step_by(10_000)
            .filter(|x| (*x < 4_000_000) || (*x > 6_000_000))
            .collect();
        let mut expected = Intervals::<u32>::new();
        expected.push(4_000_000..6_000_000);
        (genome, pos, expected)
    }

    #[test]
    fn test_find_intervals_with_low_snps_density() {
        let (genome, pos, expected) = prepare_data();
        let observed =
            find_intervals_with_low_snps_density(1.0, 10.0, &pos, genome.ginfo(), genome.gmap());
        assert!(observed.is_equal(&expected));

        let observed =
            find_intervals_with_low_snps_density(1.0, 200.0, &pos, genome.ginfo(), genome.gmap());
        let mut expected = Intervals::<u32>::new();
        expected.push(000_000..30_000_000);
        assert!(observed.is_equal(&expected));
    }
}
