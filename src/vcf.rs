use ahash::AHashSet;

use crate::{
    genome::GenomeInfo,
    genotype::{common::*, rare::*},
    indiv::Individuals,
    site::*,
};
use std::path::Path;

pub fn read_vcf(
    target_samples: &AHashSet<String>,
    gconfig: &GenomeInfo,
    vcf_path: impl AsRef<Path>,
    max_maf: f64,
    region: Option<(u32, u64, Option<u64>)>,
) -> (Sites, Individuals, GenotypeRecords) {
    use rust_htslib::bcf::{record::GenotypeAllele, IndexedReader, Read};
    // let  vcf_fn =  "/local/chib/TOPMedFull/data/org/87668/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.10b/phased/freeze.10b.chr22.pass_only.phased.bcf";
    // let vcf_fn = "./testdata/test.bcf";
    let mut bcf = IndexedReader::from_path(vcf_path).unwrap();

    match region {
        Some((rid, start, end_opt)) => {
            let chrname = &gconfig.chromnames[rid as usize];
            let rid2 = bcf.header().name2rid(chrname.as_bytes()).unwrap();
            bcf.fetch(rid2, start, end_opt).unwrap();
        }
        None => {}
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
    let individuals = Individuals::from_iter(
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
    for (_i, record_result) in bcf.records().enumerate() {
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

        let is_new_pos;
        if pos_last.is_none() {
            is_new_pos = true
        } else {
            let rid_last = rid_last.unwrap();
            let pos_last = pos_last.unwrap();
            if rid_last == chrid {
                assert!(pos_last <= pos);
            }
            if (pos_last == pos) && (rid_last == chrid) {
                is_new_pos = false;
            } else {
                is_new_pos = true;
            }
        }

        let start_allele: u8;
        if is_new_pos {
            // output information from old pos
            if pos_last.is_some() {
                let pos_last = pos_last.unwrap();
                output_rare_allele_records(
                    &mut sites,
                    &mut records,
                    &mut ab,
                    &mut allele_counts,
                    &gt,
                    pos_last,
                    ac_thres,
                );
            }

            ab.clear();
            start_allele = 0;
            ab.push(b" REF");
            for allele in alleles.iter().skip(1) {
                ab.push(b" ");
                ab.push_to_data_only(alleles[0]);
                ab.push_to_data_only(b">");
                ab.push_to_data_only(*allele);
            }
        } else {
            // Minus 1 as only non ref are added
            start_allele = ab.len() as u8 - 1;
            // skip the ref allele as it has been added
            let alleles = record.alleles();
            for allele in record.alleles().iter().skip(1) {
                ab.push(b" ");
                ab.push_to_data_only(alleles[0]);
                ab.push_to_data_only(b">");
                ab.push_to_data_only(*allele);
            }
        }
        // println!(
        // "line: {i}, pos: {pos}, start_allele: {start_allele}, ab: {:?}",
        // ab
        // );

        let fmt = record.format(b"GT");
        let bcf_fmt = fmt.inner();
        assert_eq!(bcf_fmt.type_, 1); // uint8_t
        assert_eq!(bcf_fmt.n, 2);
        assert_eq!(bcf_fmt.p_len, 2 * nsam);

        let raw_gt_slice = unsafe { std::slice::from_raw_parts(bcf_fmt.p, bcf_fmt.p_len as usize) };

        // make a vector allele idx
        let mut hap_mask = vec![];
        for sm in sample_mask.iter() {
            hap_mask.push(*sm);
            hap_mask.push(*sm);
        }

        let gt_iter = raw_gt_slice
            .iter()
            .zip(hap_mask.iter())
            // skip 'no' haplotypes
            .filter(|(_, &m)| m)
            .map(|(x, _)| {
                let mut i = match ((*x) as i32).into() {
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

        // save to ab_last in case some vcf file are normalized into biallelic lines
        // where each multi-allelic sites are encodes multiple biallelic lines
        // ab_last allows to refer alleles the previous lines
        rid_last = Some(chrid);
        pos_last = Some(pos);
    }

    // output the last record:
    if pos_last.is_some() {
        let pos_last = pos_last.unwrap();
        output_rare_allele_records(
            &mut sites,
            &mut records,
            &mut ab,
            &mut allele_counts,
            &gt,
            pos_last,
            ac_thres,
        );
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

    match region {
        Some((rid, start, end_opt)) => {
            let chrname = &gconfig.chromnames[rid as usize];
            let rid2 = bcf.header().name2rid(chrname.as_bytes()).unwrap();
            bcf.fetch(rid2, start, end_opt).unwrap();
        }
        None => {}
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
    let individuals = Individuals::from_iter(
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
    for (_i, record_result) in bcf.records().enumerate() {
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
        let _is_new_pos;
        if pos_last.is_none() {
            _is_new_pos = true
        } else {
            let rid_last = rid_last.unwrap();
            let pos_last = pos_last.unwrap();
            if rid_last == chrid {
                assert!(pos_last <= pos);
            }
            if (pos_last == pos) && (rid_last == chrid) {
                _is_new_pos = false;
            } else {
                _is_new_pos = true;
            }
        }

        let fmt = record.format(b"GT");
        let bcf_fmt = fmt.inner();
        assert_eq!(bcf_fmt.type_, 1); // uint8_t
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
            .iter()
            .zip(hap_mask.iter())
            // skip 'no' haplotypes
            .filter(|(_, &m)| m)
            .map(|(x, _)| {
                let i = match (*x as i32).into() {
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
        sites.append_bytes_to_last_allele(b">");
        sites.append_bytes_to_last_allele(alleles[1]);

        rid_last = Some(chrid);
        pos_last = Some(pos);
    }

    (sites, individuals, gm)
}

fn output_rare_allele_records(
    sites: &mut Sites,
    records: &mut Vec<GenotypeRecord>,
    ab: &mut AlleleBuffer,
    allele_counts: &mut Vec<usize>,
    gt: &Vec<u8>,
    pos: u32,
    ac_thres: usize,
) {
    // Get allele counts
    // init
    allele_counts.clear();
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

    // Set new allele encoding according the allele counts
    // new enc of rare variants  starts from 1
    let mut new_encode = 1u8;
    for (i, ac) in allele_counts.iter().enumerate() {
        assert!(i < u8::MAX as usize - 1);
        // for rare alleles
        if *ac <= ac_thres {
            ab.set_enc(i as u8, new_encode);
            new_encode += 1;
        }
        // for common alleles
        else {
            ab.set_enc(i as u8, 0);
        }
    }

    if new_encode == 1 {
        // no rare variant in this site, skip adding the variants and sites
        return;
    }

    // Remap alleleic encodings and generate rare variants records
    // all commmon alleles are encoded zeros
    // all rare allele are encoded nonzeros
    for (i, old_enc) in gt.iter().enumerate() {
        let new_enc = ab.get_enc(*old_enc as u8);
        if new_enc == 0 {
            continue;
        }
        let mut rec = GenotypeRecord::new(0);
        rec.set(pos as u32, i as u32, new_enc as u8);
        records.push(rec);
    }

    // add sites info
    sites.add_site(pos, ab);

    // if pos == 35702311 {
    //     println!(
    //         "pos={pos}, alllee={:?}, ab={:?}, ac={:?}",
    //         sites.get_site_by_position(pos),
    //         ab,
    //         allele_counts,
    //     );
    //     std::process::exit(-1);
    // }
}
