use super::{GtencodeError, Result};
use error_stack::*;

use super::Commands;
use ishare::{
    error::PathAttachment,
    genome::{Genome, GenomeInfo},
    genotype::{common::GenotypeMatrix, rare::GenotypeRecords},
    indiv::Individuals,
    site::Sites,
    utils::path::from_prefix,
    vcf::{read_vcf, read_vcf_for_genotype_matrix},
};
use rayon::prelude::*;

use rust_htslib::bcf::Read;

pub fn main_encode(args: &Commands) -> Result<()> {
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
            .change_context(GtencodeError::Input)?
            .trim()
            .split("\n")
            .for_each(|x| {
                target_samples.insert(x.to_owned());
            });
    }

    // encoding
    let ginfo = if matches!( genome_info.as_path().extension(), Some(ext) if ext == "toml") {
        GenomeInfo::from_toml_file(genome_info).change_context(GtencodeError::Input)?
    } else {
        let genome = Genome::load_from_bincode_file(genome_info.to_string_lossy().as_ref())
            .change_context(GtencodeError::Input)?;
        genome.ginfo().clone()
    };

    // divide genome into 10Mb chunks
    let mut regions = ginfo.partition_genome(parallel_chunksize_bp.map(|x| x as u32));
    // filter region with no records
    use rust_htslib::bcf::IndexedReader;
    let mut ireader = IndexedReader::from_path(vcf)
        .attach_with(|| PathAttachment::from(vcf))
        .change_context(GtencodeError::Input)?;
    let mut rec = ireader.empty_record();

    let mut nfail = 0;
    regions.retain(|r| match r.as_ref() {
        None => true,
        Some(r) => {
            let chrname = &ginfo.chromnames[r.0 as usize];
            let rid2 = match ireader.header().name2rid(chrname.as_bytes()) {
                Ok(rid2) => rid2,
                _ => {
                    // ignore regions of which chromsome name not present in the vcf header
                    return false;
                }
            };
            if ireader.fetch(rid2, r.1, r.2).is_err() {
                nfail += 1;
            }
            ireader.read(&mut rec).is_some()
        }
    });
    if nfail > 0 {
        bail!(GtencodeError::Library
            .into_report()
            .attach("region filter error"));
    }

    // construct output file names
    let gt_file = if *matrix {
        from_prefix(out, "mat").change_context(GtencodeError::Input)?
    } else {
        from_prefix(out, "rec").change_context(GtencodeError::Input)?
    };
    let sites_file = from_prefix(out, "sit").change_context(GtencodeError::Input)?;
    let ind_file = from_prefix(out, "ind").change_context(GtencodeError::Input)?;

    if *matrix {
        // parallel running
        let mut res: Vec<(Sites, Individuals, GenotypeMatrix)> = regions
            .into_par_iter()
            .map(|region| -> Result<_> {
                // (sites, individuals, GenotypeMatrix)
                let res = read_vcf_for_genotype_matrix(
                    &target_samples,
                    &ginfo,
                    vcf,
                    *threshold_maf,
                    region,
                )
                .change_context(GtencodeError::Library)?;
                if region.is_some() {
                    println!("{region:?}");
                }
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        // merge results
        let (mut sites, individuals, mut gm) = res
            .pop()
            .ok_or(GtencodeError::Library)
            .attach("Empty of results of regions")?;

        for (ss, _, rr) in res {
            sites.merge(ss);
            gm.merge(rr);
        }

        // sort
        let orders = sites
            .sort_by_position_then_allele()
            .change_context(GtencodeError::Library)?;
        let gm_ordered = gm
            .reorder_rows(&orders)
            .change_context(GtencodeError::Library)?;

        // write to files

        gm_ordered
            .into_parquet_file(&gt_file)
            .change_context(GtencodeError::Output)?;
        sites
            .into_parquet_file(&sites_file)
            .change_context(GtencodeError::Output)?;
        individuals
            .into_parquet_file(&ind_file)
            .change_context(GtencodeError::Output)?;
    } else {
        // parallel running
        let mut res: Vec<(Sites, Individuals, GenotypeRecords)> = regions
            // let mut res: Vec<_> = regions
            .into_par_iter()
            .map(|region| {
                // (sites, individuals, GenotypeRecords)
                let res = read_vcf(&target_samples, &ginfo, vcf, *threshold_maf, region)
                    .change_context(GtencodeError::Library)?;
                if region.is_some() {
                    println!("{region:?}");
                }
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        // merge results
        let (mut sites, individuals, mut records) = res
            .pop()
            .ok_or(GtencodeError::Library)
            .attach("empty of results of regions")?;

        for (ss, _, rr) in res {
            sites.merge(ss);
            records.merge(rr);
        }

        // sort
        _ = sites.sort_by_position_then_allele();
        records
            .sort_by_genome()
            .change_context(GtencodeError::Library)?;
        // write to files

        records
            .into_parquet_file(&gt_file)
            .change_context(GtencodeError::Output)?;
        sites
            .into_parquet_file(&sites_file)
            .change_context(GtencodeError::Output)?;
        individuals
            .into_parquet_file(&ind_file)
            .change_context(GtencodeError::Output)?;
    }

    // report encoding time used
    let duration = start.elapsed();
    println!("# Encoding Time : {duration:?}");
    Ok(())
}
