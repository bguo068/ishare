use super::super::Commands;
use arrow::array::{ArrayRef, UInt32Array};
use arrow::record_batch::RecordBatch;
use ishare::indiv::Individuals;
use ishare::{genotype::rare::GenotypeRecords, utils::path::from_prefix};
use itertools::{EitherOrBoth, Itertools};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use slice_group_by::GroupByMut;
use std::sync::Arc;

pub fn main_samplediff(args: &Commands) {
    if let Commands::SampleDiff { rec, pairs, out } = args {
        let mut records = GenotypeRecords::from_parquet_file(rec);
        assert!(records.is_sorted_by_genome());

        let ind_file = rec.with_extension("ind");
        let inds = Individuals::from_parquet_file(&ind_file);

        let mut res_vec = prepare_pairs(&inds, pairs);

        let nhaps = inds.m().len() * 2;
        let possible_pairs = nhaps * (nhaps - 1) / 2;
        if res_vec.len() > possible_pairs / 10 {
            // this method consolidate haplotype-level gt to sample-level gt
            // quicker for large number paris with the overhead of updating the records first
            update_res_by_consolidate_haplotypes(&mut records, &mut res_vec);
        } else {
            // no overhead but need to access four genome for each sample pair
            // quicker for smaller number of pairs
            update_res_by_compare_4_haploytpes(&records, &mut res_vec);
        }
        output(out, &res_vec, &inds)
    }
}

use std::path::PathBuf;
fn prepare_pairs(inds: &Individuals, pairs: &Option<PathBuf>) -> Vec<(u32, u32, u32)> {
    // preare pairs
    //
    let mut res_vec = Vec::<(u32, u32, u32)>::new();
    let ind_m = inds.m();
    match pairs.as_ref() {
        Some(pair_path) => std::fs::read_to_string(pair_path)
            .expect("cannot read pairs into string")
            .trim()
            .split("\n")
            .for_each(|line| {
                let mut fields = line.split(" ");
                let sample_name1 = fields.next().expect("parse 1st field");
                let sample_name2 = fields.next().expect("parse 2nd field");
                res_vec.push((ind_m[sample_name1] as u32, ind_m[sample_name2] as u32, 0));
            }),
        None => {
            for i in 0..(ind_m.len() as u32) {
                for j in 0..(ind_m.len() as u32) {
                    if i < j {
                        continue;
                    }
                    res_vec.push((i, j, 0));
                }
            }
        }
    }
    res_vec
}

fn update_res_by_compare_4_haploytpes(
    records: &GenotypeRecords,
    res_vec: &mut Vec<(u32, u32, u32)>,
) {
    res_vec.par_iter_mut().for_each(|item| {
        // for sample 1
        let it1 = {
            let g1 = item.0 * 2;
            let g2 = g1 + 1;
            records
                .iter_genome_pair_genotypes(g1, g2)
                .map(|(pos, a, b)| match (a, b) {
                    (Some(a), Some(b)) => {
                        assert_eq!(a, b, "some site is not biallelic");
                        (pos, 2)
                    }
                    (_, _) => (pos, 1),
                })
        };
        // for sample 1
        let it2 = {
            let g1 = item.1 * 2;
            let g2 = g1 + 1;
            records
                .iter_genome_pair_genotypes(g1, g2)
                .map(|(pos, a, b)| match (a, b) {
                    (Some(a), Some(b)) => {
                        assert_eq!(a, b, "some site is not biallelic");
                        (pos, 2)
                    }
                    (_, _) => (pos, 1),
                })
        };

        let mut discord = 0;
        it1.merge_join_by(it2, |a, b| a.0.cmp(&b.0))
            .for_each(|x| match x {
                EitherOrBoth::Both(a, b) => {
                    if a.1 != b.1 {
                        discord += 1;
                    }
                }
                _ => {
                    discord += 1;
                }
            });
        item.2 = discord;
    });
}
fn update_res_by_consolidate_haplotypes(
    records: &mut GenotypeRecords,
    res_vec: &mut Vec<(u32, u32, u32)>,
) {
    // Work-around:
    // considate each two genomes of of the sample into a single genome
    // ensure all sites are biallelic
    // the allele byte will be encode as count of rare allele per site instead of allele id
    let v = records.records_mut();
    // . resort
    v.sort_by_key(|r| {
        let (pos, genome, allele) = r.get();
        (genome - genome % 2, pos, allele)
    });
    // . combine genotype and mark reduant record for deletion
    v.as_mut_slice()
        .linear_group_by_key_mut(|r| {
            let (pos, genome, _allele) = r.get();
            (genome - genome % 2, pos)
        })
        .for_each(|grp| match grp.len() {
            1 => {
                let (pos, mut genome, _allele) = grp[0].get();
                genome = genome - genome % 2;
                grp[0].set(pos, genome, 1);
            }
            2 => {
                let (pos, mut genome, allele1) = grp[0].get();
                let (_, _, allele2) = grp[1].get();
                // check if biallelic sites: if there two rare allele, they should the same
                assert_eq!(allele1, allele2, "it looks like some site is not biallelic");
                genome = genome - genome % 2;
                grp[0].set(pos, genome, 2);
                grp[1].set_sentinel();
            }
            _ => {
                panic!("more than two allele at a give sample and site");
            }
        });
    // clean up empty values
    v.retain(|r| !r.is_sentinel());

    res_vec.par_iter_mut().for_each(|item| {
        let g1 = item.0 * 2;
        let g2 = item.1 * 2;
        let mut discord = 0;
        records
            .iter_genome_pair_genotypes(g1, g2)
            .for_each(|(_pos, a, b)| match (a, b) {
                (Some(a), Some(b)) if a == b => {}
                (_, _) => {
                    discord += 1;
                }
            });
        item.2 = discord;
    });
}

fn output(out: &Option<String>, res_vec: &[(u32, u32, u32)], inds: &Individuals) {
    // write parquet files
    match out.as_ref() {
        Some(output) => {
            println!("WARN: output option is specified, results are not printed on the screen");
            let out_path = from_prefix(&output, "pq").unwrap();
            let sample1 = UInt32Array::from(res_vec.iter().map(|item| item.0).collect_vec());
            let sample2 = UInt32Array::from(res_vec.iter().map(|item| item.1).collect_vec());
            let discord = UInt32Array::from(res_vec.iter().map(|item| item.2).collect_vec());

            let batch = RecordBatch::try_from_iter(vec![
                ("sample_id1", Arc::new(sample1) as ArrayRef),
                ("sample_id2", Arc::new(sample2) as ArrayRef),
                ("discordance", Arc::new(discord) as ArrayRef),
            ])
            .unwrap();

            let file = std::fs::File::create(&out_path).unwrap();
            let props = WriterProperties::builder()
                .set_compression(Compression::SNAPPY)
                .build();
            let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();

            writer.write(&batch).expect("writing batch");

            writer.close().unwrap();
        }
        None => {
            for (s1, s2, discordance) in res_vec {
                print!(
                    "{}\t{}\t{}\n",
                    inds.v()[*s1 as usize].as_str(),
                    inds.v()[*s2 as usize].as_str(),
                    *discordance
                );
            }
        }
    }
}
