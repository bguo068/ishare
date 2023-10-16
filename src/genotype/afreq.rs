use ahash::{HashMap, HashMapExt};
use rust_htslib::{
    bcf::{self, record::GenotypeAllele, Read},
    htslib::bcf_is_snp,
};

use crate::genome::GenomeInfo;

/// similar to get_afreq_from_vcf but use genome-wide coordinates
///
/// return a vector of tuple with element being a tuple of gw_pos and afreq
pub fn get_afreq_from_vcf_genome_wide(
    vcf_files: &Vec<String>,
    ginfo: &GenomeInfo,
) -> Vec<(u32, f32)> {
    let res_map = get_afreq_from_vcf(vcf_files);
    let mut res_vec = Vec::<(u32, f32)>::new();
    println!("number chromosome: {}", res_vec.len());
    for (chrname, v) in res_map.iter() {
        let chrid = ginfo.idx[chrname];
        for (pos, afreq) in v.iter() {
            let gw_pos = ginfo.to_gw_pos(chrid, *pos);
            res_vec.push((gw_pos, *afreq));
        }
    }
    res_vec.sort_by_key(|tuple| tuple.0);
    res_vec
}

/// get allele frequency from a list of vcf files
///
/// Only biallic SNPs are considerred; others are silently ignored
///
/// the returned results is a hashmap, with key being chromosome name,
/// and the value for each chromosme is a vector of 2-tuple (chromosomal bp position and the allele frequency)
pub fn get_afreq_from_vcf(vcf_files: &Vec<String>) -> HashMap<String, Vec<(u32, f32)>> {
    let mut res = HashMap::<String, Vec<(u32, f32)>>::new();

    for vcf_file in vcf_files.iter() {
        let mut last_chr_id = u32::MAX;
        let mut last_v = &mut vec![];
        let mut reader =
            bcf::Reader::from_path(vcf_file).expect(&format!("cannot open vcf file {}", vcf_file));
        let header = reader.header().clone();
        let mut rec = reader.empty_record();
        while let Some(_) = reader.read(&mut rec) {
            // only consider biallelic SNPs
            if !((rec.allele_count() == 2) && unsafe { bcf_is_snp(rec.inner) != 0 }) {
                continue;
            }
            let mut ref_alt_cnt = [0u32; 2];
            for gt_sample in rec
                .format(b"GT")
                .integer()
                .expect("failure to get genotype")
                .iter()
            {
                for gta in *gt_sample {
                    let a: GenotypeAllele = (*gta).into();
                    match a.index() {
                        Some(i) => ref_alt_cnt[i as usize] += 1,
                        None => {}
                    }
                }
            }

            let total = ref_alt_cnt[0] + ref_alt_cnt[1];
            if total == 0 {
                continue;
            }
            let afreq = ref_alt_cnt[1] as f32 / (total as f32);

            let pos = rec.pos() as u32;
            let chrid = rec.rid().unwrap();
            if chrid == last_chr_id {
                last_v.push((pos, afreq));
            } else {
                let last_chr_name = std::str::from_utf8(header.rid2name(chrid).unwrap()).unwrap();
                println!("chrname: {}", last_chr_name);
                res.entry(last_chr_name.to_owned())
                    .and_modify(|v| v.push((pos, afreq)))
                    .or_insert(vec![(pos, afreq)]);
                last_chr_id = chrid;
                last_v = res.get_mut(last_chr_name).unwrap();
            }
        }
    }
    res
}
