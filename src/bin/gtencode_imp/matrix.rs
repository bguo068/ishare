use ishare::{genotype::common::GenotypeMatrix, site::Sites};
use itertools::Itertools;

use super::super::Commands;

pub fn main_matrix(args: &Commands) {
    if let Commands::Matrix {
        mat,
        genomes,
        positions,
    } = args
    {
        let gm = GenotypeMatrix::from_parquet_file(mat);
        let sit_file = mat.clone().with_extension("sit");
        let sites = Sites::from_parquet_file(sit_file);

        let genomes = match genomes {
            Some(v) if v.len() == 0 => Vec::new(),
            Some(v) => v.clone(),
            None => Vec::new(),
        };
        let pos_idx = match positions {
            Some(v) if v.len() == 0 => Vec::new(),
            Some(v) => {
                let mut u = Vec::<u32>::new();
                for pos in v {
                    let (s, e) = sites.get_idx_by_position(*pos);
                    if s == e {
                        eprintln!("position {} not found in the matrix", pos);
                        panic!("position not found");
                    }
                    u.extend((s..e).map(|x| x as u32));
                }

                u
            }
            None => Vec::new(),
        };

        let nrows = gm.nrows();
        // let ncols: usize = gm.get_ncol();

        use bstr::ByteSlice;
        match (pos_idx.len() == 0, genomes.len() == 0) {
            (true, true) => {
                for row in 0..nrows {
                    let (pos, bytes) = sites.get_site_by_idx(row as usize);
                    print!("pos={:10}, allele={:10.10}\t", pos, bytes.as_bstr());
                    let s = gm
                        .get_row(row)
                        .iter()
                        .map(|x| if *x { "1" } else { "0" })
                        .join("");
                    println!("{s}");
                }
            }
            (true, false) => {
                for row in 0..nrows {
                    let (pos, bytes) = sites.get_site_by_idx(row as usize);
                    print!("pos={:10}, allele={:10.10}\t", pos, bytes.as_bstr());
                    let s = genomes
                        .iter()
                        .map(|x| {
                            let hit = gm.get_at(row as usize, *x as usize);
                            if hit {
                                "1"
                            } else {
                                "0"
                            }
                        })
                        .join("");
                    println!("{s}");
                }
            }
            (false, true) => {
                for row in pos_idx {
                    let (pos, bytes) = sites.get_site_by_idx(row as usize);
                    print!("pos={:10}, allele={:10.10}\t", pos, bytes.as_bstr());
                    let s = gm
                        .get_row(row as usize)
                        .iter()
                        .map(|x| if *x { "1" } else { "0" })
                        .join("");
                    println!("{s}");
                }
            }
            (false, false) => {
                for row in pos_idx {
                    let (pos, bytes) = sites.get_site_by_idx(row as usize);
                    print!("pos={:10}, allele={:10.10}\t", pos, bytes.as_bstr());
                    let s = genomes
                        .iter()
                        .map(|x| {
                            let hit = gm.get_at(row as usize, *x as usize);
                            if hit {
                                "1"
                            } else {
                                "0"
                            }
                        })
                        .join("");
                    println!("{s}");
                }
            }
        }
    }
}
