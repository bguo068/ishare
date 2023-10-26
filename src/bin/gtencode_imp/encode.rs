use super::super::Commands;
use ishare::{
    genome::GenomeInfo,
    utils::path::from_prefix,
    vcf::{read_vcf, read_vcf_for_genotype_matrix},
};
use rayon::prelude::*;
use rust_htslib::bcf::Read;

pub fn main_encode(args: &Commands) {
    // unpack cli args
    let (vcf, sample_lst, genome_info, out, parallel_chunksize_bp, matrix, threshold_maf) =
        if let Commands::Encode {
            vcf,
            samples_lst,
            genome_info,
            out,
            parallel_chunksize_bp,
            matrix,
            threshold_maf,
        } = args
        {
            (
                vcf,
                samples_lst,
                genome_info,
                out,
                parallel_chunksize_bp,
                matrix,
                threshold_maf,
            )
        } else {
            panic!("wrong type")
        };

    use std::time::Instant;
    let start = Instant::now();
    println!("# Encoding genotypes ...");

    // read sample list
    use ahash::AHashSet;
    let mut target_samples = AHashSet::new();
    if let Some(sample_lst) = sample_lst.as_ref() {
        std::fs::read_to_string(sample_lst)
            .expect("can not open sample list file")
            .trim()
            .split("\n")
            .for_each(|x| {
                target_samples.insert(x.to_owned());
            });
    }

    // encoding
    let ginfo = GenomeInfo::from_toml_file(genome_info);

    // divide genome into 10Mb chunks
    let mut regions = ginfo.partition_genome(parallel_chunksize_bp.map(|x| x as u32));
    // filter region with no records
    use rust_htslib::bcf::IndexedReader;
    let mut ireader = IndexedReader::from_path(vcf).unwrap();
    let mut rec = ireader.empty_record();
    regions.retain(|r| match r.as_ref() {
        None => true,
        Some(r) => {
            let chrname = &ginfo.chromnames[r.0 as usize];
            let rid2 = ireader.header().name2rid(chrname.as_bytes()).unwrap();
            ireader.fetch(rid2, r.1, r.2).unwrap();
            match ireader.read(&mut rec) {
                None => false,
                Some(_) => true,
            }
        }
    });

    // construct output file names
    let gt_file = if *matrix {
        from_prefix(out, "mat").unwrap()
    } else {
        from_prefix(out, "rec").unwrap()
    };
    let sites_file = from_prefix(out, "sit").unwrap();
    let ind_file = from_prefix(out, "ind").unwrap();

    if *matrix {
        // parallel running
        let mut res: Vec<_> = regions
            .into_par_iter()
            .map(|region| {
                // (sites, individuals, GenotypeMatrix)
                let res = read_vcf_for_genotype_matrix(
                    &target_samples,
                    &ginfo,
                    vcf,
                    *threshold_maf,
                    region,
                );
                if region.is_some() {
                    println!("{:?}", region);
                }
                res
            })
            .collect();

        // merge results
        let (mut sites, individuals, mut gm) = res.pop().unwrap();

        for (ss, _, rr) in res {
            sites.merge(ss);
            gm.merge(rr);
        }

        // sort
        let orders = sites.sort_by_position_then_allele();
        let gm_ordered = gm.reorder_rows(&orders);

        // write to files

        gm_ordered.into_parquet_file(&gt_file);
        sites.into_parquet_file(&sites_file);
        individuals.into_parquet_file(&ind_file);
    } else {
        // parallel running
        let mut res: Vec<_> = regions
            .into_par_iter()
            .map(|region| {
                // (sites, individuals, GenotypeRecords)
                let res = read_vcf(&target_samples, &ginfo, vcf, *threshold_maf, region);
                if region.is_some() {
                    println!("{:?}", region);
                }
                res
            })
            .collect();

        // merge results
        let (mut sites, individuals, mut records) = res.pop().unwrap();

        for (ss, _, rr) in res {
            sites.merge(ss);
            records.merge(rr);
        }

        // sort
        _ = sites.sort_by_position_then_allele();
        records.sort_by_genome();
        // write to files

        records.into_parquet_file(&gt_file);
        sites.into_parquet_file(&sites_file);
        individuals.into_parquet_file(&ind_file);
    }

    // report encoding time used
    let duration = start.elapsed();
    println!("# Encoding Time : {:?}", duration);
}
