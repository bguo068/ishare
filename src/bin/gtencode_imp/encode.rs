use super::super::Commands;
use ishare::{
    genome::GenomeInfo,
    vcf::{read_vcf, read_vcf_for_genotype_matrix},
};
use rayon::prelude::*;
use std::path::{Path, PathBuf};

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
    let regions = ginfo.partition_genome(parallel_chunksize_bp.map(|x| x as u32));

    // construct output file names
    let dir = out.parent().unwrap().to_str().unwrap();
    if !Path::new(dir).exists() {
        std::fs::create_dir_all(dir).unwrap();
    }
    let mut prefix = out.file_name().unwrap().to_owned().into_string().unwrap();
    if prefix.ends_with(".rec") {
        prefix = prefix.strip_suffix(".rec").unwrap().to_string();
    }
    if prefix.ends_with(".mat") {
        prefix = prefix.strip_suffix(".mat").unwrap().to_string();
    }
    let mut gt_file = PathBuf::from(dir);
    if *matrix {
        gt_file.push(PathBuf::from(format!("{prefix}.mat")));
    } else {
        gt_file.push(PathBuf::from(format!("{prefix}.rec")));
    }
    let mut sites_file = PathBuf::from(dir);
    sites_file.push(PathBuf::from(format!("{prefix}.sit")));
    let mut ind_file = PathBuf::from(dir);
    ind_file.push(PathBuf::from(format!("{prefix}.ind")));

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
