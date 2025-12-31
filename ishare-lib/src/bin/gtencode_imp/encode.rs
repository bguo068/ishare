use std::backtrace::Backtrace;

use super::super::Commands;
use ishare::{
    genome::{Genome, GenomeInfo},
    genotype::{common::GenotypeMatrix, rare::GenotypeRecords},
    indiv::Individuals,
    site::Sites,
    utils::path::from_prefix,
    vcf::{read_vcf, read_vcf_for_genotype_matrix},
};
use rayon::prelude::*;

use rust_htslib::bcf::Read;

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    // outside
    #[snafu(transparent)]
    Genome {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genome::Error,
    },
    #[snafu(transparent)]
    UtilsPath {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::utils::path::Error,
    },
    #[snafu(transparent)]
    Vcf {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::vcf::Error,
    },

    #[snafu(transparent)]
    Sites {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::site::Error,
    },
    #[snafu(transparent)]
    GenotypeCommon {
        // non leaf
        #[snafu(backtrace)]
        #[snafu(source(from(ishare::genotype::common::Error, Box::new)))]
        source: Box<ishare::genotype::common::Error>,
    },
    #[snafu(transparent)]
    GenotypeRare {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genotype::rare::Error,
    },
    #[snafu(transparent)]
    Individual {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::indiv::Error,
    },

    // this module
    StdIo {
        // leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    Hts {
        // leaf
        #[snafu(source(from(rust_htslib::errors::Error, Box::new)))]
        source: Box<rust_htslib::errors::Error>,
        backtrace: Box<Option<Backtrace>>,
    },
    RegionFilter {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
    EmptyVec {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
}
type Result<T> = std::result::Result<T, Error>;

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
            .context(StdIoSnafu {})?
            .trim()
            .split("\n")
            .for_each(|x| {
                target_samples.insert(x.to_owned());
            });
    }

    // encoding
    let ginfo = if matches!( genome_info.as_path().extension(), Some(ext) if ext == ".toml") {
        GenomeInfo::from_toml_file(genome_info)?
    } else {
        let genome = Genome::load_from_bincode_file(genome_info.to_string_lossy().as_ref())?;
        genome.ginfo().clone()
    };

    // divide genome into 10Mb chunks
    let mut regions = ginfo.partition_genome(parallel_chunksize_bp.map(|x| x as u32));
    // filter region with no records
    use rust_htslib::bcf::IndexedReader;
    let mut ireader = IndexedReader::from_path(vcf).context(HtsSnafu {})?;
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
        RegionFilterSnafu {}.fail()?;
    }

    // construct output file names
    let gt_file = if *matrix {
        from_prefix(out, "mat")?
    } else {
        from_prefix(out, "rec")?
    };
    let sites_file = from_prefix(out, "sit")?;
    let ind_file = from_prefix(out, "ind")?;

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
                )?;
                if region.is_some() {
                    println!("{region:?}");
                }
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        // merge results
        let (mut sites, individuals, mut gm) = res.pop().context(EmptyVecSnafu)?;

        for (ss, _, rr) in res {
            sites.merge(ss);
            gm.merge(rr);
        }

        // sort
        let orders = sites.sort_by_position_then_allele()?;
        let gm_ordered = gm.reorder_rows(&orders)?;

        // write to files

        gm_ordered.into_parquet_file(&gt_file)?;
        sites.into_parquet_file(&sites_file)?;
        individuals.into_parquet_file(&ind_file)?;
    } else {
        // parallel running
        let mut res: Vec<(Sites, Individuals, GenotypeRecords)> = regions
            // let mut res: Vec<_> = regions
            .into_par_iter()
            .map(|region| {
                // (sites, individuals, GenotypeRecords)
                let res = read_vcf(&target_samples, &ginfo, vcf, *threshold_maf, region)?;
                if region.is_some() {
                    println!("{region:?}");
                }
                Ok(res)
            })
            .collect::<Result<Vec<_>>>()?;

        // merge results
        let (mut sites, individuals, mut records) = res.pop().context(EmptyVecSnafu {})?;

        for (ss, _, rr) in res {
            sites.merge(ss);
            records.merge(rr);
        }

        // sort
        _ = sites.sort_by_position_then_allele();
        records.sort_by_genome()?;
        // write to files

        records.into_parquet_file(&gt_file)?;
        sites.into_parquet_file(&sites_file)?;
        individuals.into_parquet_file(&ind_file)?;
    }

    // report encoding time used
    let duration = start.elapsed();
    println!("# Encoding Time : {duration:?}");
    Ok(())
}
