use ishare::{
    genome::GenomeInfo, genotype::rare::GenotypeRecords, indiv::Individuals, site::Sites,
    utils::path::from_prefix,
};
use parquet::data_type::AsBytes;
use rust_htslib::bcf::{self, record::GenotypeAllele};
use slice_group_by::*;

use itertools::EitherOrBoth::*;
use itertools::Itertools;

use super::args::Commands;

pub fn main_export(args: &Commands) {
    if let Commands::Export {
        rec,
        output_fmt,
        genome_info,
        out_prefix,
        keep_sites_with_multi_common_alleles,
    } = args
    {
        let mut records = GenotypeRecords::from_parquet_file(&from_prefix(rec, "rec").unwrap());
        let sites = Sites::from_parquet_file(&from_prefix(rec, ".sit").unwrap());
        let inds = Individuals::from_parquet_file(&from_prefix(rec, "ind").unwrap());
        let ginfo = GenomeInfo::from_toml_file(genome_info);

        let mut compressed = true;
        let mut format = bcf::Format::Bcf;
        match output_fmt.as_str() {
            "b" => {}
            "u" => compressed = false,
            "z" => format = bcf::Format::Vcf,
            "v" => {
                compressed = false;
                format = bcf::Format::Vcf
            }
            _ => {
                eprint!("unsupported output file format");
                std::process::exit(-1);
            }
        }

        // prepare bcf header
        let mut header = bcf::Header::new();
        //// 1. add contig lines
        let mut line = Vec::<u8>::new();
        use std::io::Write;
        ginfo
            .chromnames
            .as_slice()
            .iter()
            .zip(ginfo.chromsize.as_slice().iter())
            .for_each(|(chrname, chrsize)| {
                line.clear();
                write!(line, "##contig=<ID={},length={}>", chrname, chrsize).unwrap();
                header.push_record(line.as_bytes());
            });
        //// 2. add FORMAT/GT
        let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
        header.push_record(header_gt_line.as_bytes());
        //// 3. add samples
        inds.v().iter().for_each(|sample| {
            header.push_sample(sample.as_bytes());
        });

        let mut writer = match out_prefix {
            Some(out) => bcf::Writer::from_path(out, &header, !compressed, format).unwrap(),
            None => bcf::Writer::from_stdout(&header, !compressed, format).unwrap(),
        };

        let mut bcf_record = writer.empty_record();

        // sort by postion
        records.sort_by_position();

        let mut allele = Vec::new();
        let mut is_rare = Vec::new();
        let mut gt = Vec::new();
        let nsam = writer.header().sample_count();
        // write records
        records
            .records()
            .linear_group_by_key(|r| r.get_position())
            // use merge join by to efficiently get the position's index in Sites
            .merge_join_by(sites.get_gw_pos_slice().iter().enumerate(), |a, b| {
                a[0].get_position().cmp(b.1)
            })
            // records all this site
            .for_each(|mrg| {
                match mrg {
                    Left(_) => panic!("unreachable"),
                    Right(_) => panic!("unreachable"),
                    Both(recs, (ipos, _)) => {
                        let gw_pos = recs[0].get_position();
                        let (_chrid, chrname, pos) = ginfo.to_chr_pos(gw_pos);
                        let rid = writer.header().name2rid(chrname.as_bytes()).unwrap();

                        allele.clear();
                        gt.clear();
                        is_rare.clear();

                        let allele_bytes = sites.get_alleles_by_idx(ipos);
                        for a in allele_bytes.split(|byte| *byte != b' ') {
                            allele.push(a);
                        }

                        is_rare.resize(allele.len(), false);
                        recs.iter().for_each(|r| {
                            let i = r.get_allele() as usize;
                            if is_rare[i] == false {
                                is_rare[i] = true;
                            }
                        });
                        let mut common_allele = u8::MAX;
                        let num_common_alleles: usize = is_rare
                            .iter()
                            .enumerate()
                            .map(|(i, isr)| match isr {
                                false => {
                                    // get the first common allele (if there are multiple common alleles)
                                    if common_allele == u8::MAX {
                                        common_allele = i as u8;
                                    }
                                    i
                                }
                                true => 0,
                            })
                            .sum();
                        if num_common_alleles > 1 {
                            if *keep_sites_with_multi_common_alleles {
                                eprint!("use the first common allele to represent all common allele for a multi-common-allele site!");
                            } else {
                                eprint!("Skiping a site as it has multiple common alleles!");
                                return;
                            }
                        }

                        bcf_record.set_rid(Some(rid));
                        bcf_record.set_pos(pos as i64);
                        bcf_record.set_alleles(allele.as_slice()).unwrap();
                        // figure out what alleles are rare
                        // convert to genotype vector
                        recs.iter()
                            .merge_join_by(0..(2 * nsam), |r, igenome| r.get_genome().cmp(igenome))
                            .for_each(|mrg2| {
                                let gta = match mrg2 {
                                    Left(_) => panic!("unreachable"),
                                    Right(ignome) => match ignome % 2 == 0 {
                                        true => GenotypeAllele::Unphased(common_allele as i32),
                                        false => GenotypeAllele::Phased(common_allele as i32),
                                    },
                                    Both(r, ignome) => match ignome % 2 == 0 {
                                        true => GenotypeAllele::Unphased(r.get_allele() as i32),
                                        false => GenotypeAllele::Phased(r.get_allele() as i32),
                                    },
                                };
                                gt.push(gta);
                            });
                        bcf_record.push_genotypes(&gt[..]).unwrap();

                        writer.write(&bcf_record).unwrap();
                    }
                };
            });
    }
}
