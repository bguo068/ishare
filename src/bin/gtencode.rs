use ahash::{AHashMap, HashSet, HashSetExt};
use clap::{Parser, Subcommand};
use ishare::{
    genome::GenomeInfo,
    genotype::{common::GenotypeMatrix, rare::GenotypeRecords},
    indiv::Individuals,
    io::IntoParquet,
    share::mat::NamedMatrix,
    site::Sites,
    vcf::{read_vcf, read_vcf_for_genotype_matrix},
};
use itertools::Itertools;
use rayon::prelude::*;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name="gtencode", author, version, about, long_about=None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Encode genotype data from VCF file to tables or matrix (-m)
    Encode {
        /// Path to VCF (input)
        vcf: PathBuf,
        /// Path to genome info toml file (input)
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to Genotype record table (output)
        #[arg(short, long)]
        out: PathBuf,
        /// Set chunks size for parallelized parsing, if not set, the whole genome is treated as a single chunk
        #[arg(short = 'c', long)]
        parallel_chunksize_bp: Option<usize>,
        /// If set, encode genotype as matrix (not efficient for rare variants);
        /// otherwise encode genotype as table (efficient for rare variants)
        #[arg(short, long, default_value_t = false)]
        matrix: bool,
        /// For table encoding only keep variants with maf no larger than the threshold;
        /// for matrix encoding only keep variants with maf no less than the threshold;
        #[arg(short = 'T', long, default_value_t = 0.001)]
        threshold_maf: f64,
    },

    /// View encoded table records
    Records {
        /// Path to Genotype record table (input)
        rec: PathBuf,
        /// optional genome id, if not specified, show all genomes
        #[arg(short, long)]
        genome: Option<u32>,
        /// optional genome-wide position, if not specified, show all positions
        #[arg(short, long)]
        pos: Option<u32>,
    },
    /// View encoded matrix contents
    Matrix {
        mat: PathBuf,
        /// optional genome id, if not specified, show all genomes
        #[arg(short, long)]
        genomes: Option<Vec<u32>>,
        /// optional genome-wide position, if not specified, show all positions
        #[arg(short, long)]
        positions: Option<Vec<u32>>,
    },
    /// View encoded sites
    Sites {
        /// Path to Genotype record table (input)
        sit: PathBuf,
        /// optional genome-wide pos, if not specified, show all positions
        #[arg(short, long)]
        pos: Option<u32>,
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
    },
    /// View individuals
    Samples {
        /// Path to Individal file (input)
        ind: PathBuf,
        /// optional sample name if not species show all samples.
        #[arg(short, long)]
        sample: Option<String>,
        /// optional sample idx if not species show all samples.
        #[arg(short = 'I', long)]
        idx_sample: Option<usize>,
        #[arg(short = 'i', long)]
        idx_genome: Option<usize>,
    },
    /// View genotype sharing between two genome
    Share {
        /// Path to Genotype record table (input)
        rec: PathBuf,
        /// genome 1
        a: u32,
        /// genome 2
        b: u32,
    },

    /// Calculate pairwise similarity via Jaccard index
    Jaccard {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of jaccard value of a genome pair to be printed
        #[arg(short = 'j', long)]
        min_jaccard: Option<f64>,
        /// optional low threshold of total number of sites with rare alleles for a genome pair
        #[arg(short = 't', long)]
        min_total: Option<u32>,
        /// optional low threshold of total number of sites where a genome pair shares a rare allele
        #[arg(short = 's', long)]
        min_shared: Option<u32>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },

    /// Calculate pairwise similarity via Cosine
    // See definition in https://en.wikipedia.org/wiki/Cosine_similarity
    Cosine {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of consine similarity for a genome pair
        #[arg(short = 'c', long)]
        min_cosine: Option<f64>,
        /// optional low threshold of magnitude (denominator) for a genome pair
        #[arg(short = 'm', long)]
        min_magnitude: Option<f64>,
        /// optional low threshold of dot product (numerator) for a genome pair
        #[arg(short = 'p', long)]
        min_dot_prod: Option<i32>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
    /// Calculate pairwise similarity via GRM (GCTA formula)
    // See defintion in Eqn3 of GCTA paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/
    Grm {
        /// Path to genotype record table (input)
        rec: PathBuf,
        /// optional subsets of genome, if not set, use all genomes
        #[arg(short, long, group = "genome_selection")]
        genomes: Option<Vec<u32>>,
        /// optional paths to one or two genome lists (each row is a genome ids).
        /// If one list is provided, calculate within-list sharing.
        /// If two lists are provided, calculate inter-list sharing.
        /// The program will refuse to run if two overlapping lists are provided.
        #[arg(short = 'l', long, group = "genome_selection")]
        lists: Vec<PathBuf>,
        /// optional low threshold of grm relationship for a genome pair
        #[arg(short = 'r', long)]
        min_grm_related: Option<f64>,
        /// path to output file '*.jac'. If specified, results will not be printed on the screen
        #[arg(short = 'o', long)]
        output: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Encode {
            vcf,
            out,
            genome_info,
            parallel_chunksize_bp,
            matrix,
            threshold_maf,
        }) => {
            use std::time::Instant;
            let start = Instant::now();
            println!("# Encoding genotypes ...");

            // encoding
            let ginfo = GenomeInfo::from_toml_file(genome_info);

            // divide genome into 10Mb chunks
            let mut regions = Vec::new();
            match *parallel_chunksize_bp {
                Some(parallel_chunksize_bp) => {
                    for (chrid, chrlen) in ginfo.chromsize.iter().enumerate() {
                        let chrid = chrid as u32;
                        let chrlen = *chrlen;
                        let mut s = 0u32;
                        while s < chrlen {
                            let mut e = s + parallel_chunksize_bp as u32;
                            if e > chrlen {
                                e = chrlen;
                            }
                            regions.push(Some((chrid as u32, s as u64, Some(e as u64))));
                            s = e + 1;
                        }
                    }
                }
                None => {
                    regions.push(None);
                }
            };

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
                        let res = read_vcf_for_genotype_matrix(&ginfo, vcf, *threshold_maf, region);
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
                        let res = read_vcf(&ginfo, vcf, *threshold_maf, region);
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
        Some(Commands::Records { rec, genome, pos }) => {
            // read records to file
            use std::time::Instant;
            let start = Instant::now();
            println!("# Loading genotype records ...");
            let records = GenotypeRecords::from_parquet_file(rec);
            let duration = start.elapsed();
            println!("# Loading Time : {:?}", duration);
            for r in records.records() {
                let mut to_print = true;
                match pos {
                    Some(p) if r.get_position() != *p => {
                        to_print = false;
                    }
                    _ => {}
                };
                match genome {
                    Some(g) if r.get_genome() != *g => {
                        to_print = false;
                    }
                    _ => {}
                };
                if to_print {
                    println!("{:?}", r)
                }
            }
        }
        Some(Commands::Matrix {
            mat,
            genomes,
            positions,
        }) => {
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
        Some(Commands::Sites {
            sit,
            pos,
            genome_info,
        }) => {
            let ginfo = GenomeInfo::from_toml_file(genome_info);

            let sites = Sites::from_parquet_file(sit);
            for i in 0..sites.len() {
                let (p, alleles) = sites.get_site_by_idx(i);
                match pos {
                    Some(pos) if *pos != p => {
                        continue;
                    }
                    _ => {}
                }
                let alleles = std::str::from_utf8(alleles).unwrap();
                let (chrid, chrname, chrpos) = ginfo.to_chr_pos(p);
                println!(
                    "gwpos={p}, chrid={chrid}, chrname={chrname}, pos={chrpos} alleles = {}",
                    alleles
                );
            }
        }
        Some(Commands::Samples {
            ind,
            sample,
            idx_sample,
            idx_genome,
        }) => {
            let inds = Individuals::from_parquet_file(ind);

            let print_sample = |i: usize, s: &str| {
                println!(
                    "idx={i}, genome1={}, genome2={}, name={s}",
                    i * 2,
                    i * 2 + 1
                );
            };

            #[derive(Copy, Clone)]
            enum Selection {
                NONE,
                ONE(usize),
                ALL,
            }

            let mut sel = match sample {
                Some(s) => Selection::ONE(inds.m()[s]),
                None => Selection::ALL,
            };
            sel = match (idx_sample, sel) {
                (Some(i), Selection::ONE(ii)) if *i == ii => sel,
                (Some(i), Selection::ONE(ii)) if *i != ii => Selection::NONE,
                (Some(_i), Selection::NONE) => Selection::NONE,
                (Some(i), Selection::ALL) => Selection::ONE(*i),
                _ => sel,
            };
            sel = match (idx_genome, sel) {
                (Some(i), Selection::ONE(ii)) if (*i / 2) == ii => sel,
                (Some(i), Selection::ONE(ii)) if (*i / 2) != ii => Selection::NONE,
                (Some(_i), Selection::NONE) => Selection::NONE,
                (Some(i), Selection::ALL) => Selection::ONE(*i / 2),
                _ => sel,
            };
            match sel {
                Selection::NONE => {}
                Selection::ALL => {
                    for (i, s) in inds.v().iter().enumerate() {
                        print_sample(i, s);
                    }
                }
                Selection::ONE(i) => {
                    let s = &inds.v()[i];
                    print_sample(i, s);
                }
            }
        }
        Some(Commands::Share { rec, a, b }) => {
            let records = GenotypeRecords::from_parquet_file(rec);
            for (pos, gt1, gt2) in records.iter_genome_pair_genotypes(*a, *b) {
                println!("pos={}, allele_a={:?}, allelle_b={:?}", pos, gt1, gt2);
            }
        }
        Some(Commands::Jaccard {
            rec,
            genomes,
            lists,
            min_jaccard,
            min_total,
            min_shared,
            output,
        }) => {
            let min_jaccard = match min_jaccard {
                Some(x) => *x,
                None => -1.0f64,
            };
            let min_total = match min_total {
                Some(x) => *x,
                None => 0u32,
            };
            let min_shared = match min_shared {
                Some(x) => *x,
                None => 0u32,
            };

            let records = GenotypeRecords::from_parquet_file(rec);
            assert!(records.is_sorted_by_genome());

            let ind_file = rec.with_extension("ind");
            let _inds = Individuals::from_parquet_file(&ind_file);

            let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);

            // run in parallel and collect row results
            let res: Vec<(u32, u32, u32, u32)> = pairs
                .into_par_iter()
                .map(|(genome1, genome2)| {
                    let mut total: u32 = 0;
                    let mut shared: u32 = 0;
                    for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                        match (a, b) {
                            (Some(a), Some(b)) if a == b => {
                                shared += 1;
                                total += 1
                            }
                            (None, None) => {}
                            (_, _) => total += 1,
                        }
                    }

                    let mut out = Some((genome1, genome2, total, shared));
                    if (total < min_total)
                        || (shared < min_shared)
                        || ((shared as f64) / (total as f64) < min_jaccard)
                    {
                        out = None;
                    }

                    out
                })
                .flatten()
                .collect();

            // for identicial sets, elements correponsing to lower matrix is not
            // calculated use this as an indicator to fill the the low part of
            // the matrix when updating the full jaccard matrix
            let identifical = row_genomes == col_genomes;

            let mut resmat = NamedMatrix::new(row_genomes, col_genomes);

            for (g1, g2, total, shared) in res {
                let jaccard = ((shared as f32) / (total as f32) * 1e6) as u32;

                if output.is_none() {
                    println!(
                        "g1={g1}, g2={g2}, total={total}, shared={shared}, jaccard*10^6={}",
                        jaccard
                    );
                }

                // update the matrix
                resmat.set_by_names(g1, g2, jaccard);
                if identifical {
                    resmat.set_by_names(g2, g1, jaccard);
                }
            }

            // write matrix to files
            match output.as_ref() {
                Some(output) => {
                    let dir = output.parent().unwrap().to_str().unwrap();
                    let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                    if !filename.ends_with(".jac") {
                        filename.push_str(".jac");
                    }
                    let mut p = PathBuf::from(dir);
                    p.push(filename);

                    // let resmat0 = resmat.clone();
                    println!("WARN: output option is specified, results are not printed on the screen, check file {:?}", p);
                    // println!("\n writing...");
                    resmat.into_parquet(&p)
                }
                None => {}
            }
        }
        Some(Commands::Cosine {
            rec,
            genomes,
            lists,
            min_cosine,
            min_magnitude,
            min_dot_prod,
            output,
        }) => {
            let min_cosine = match min_cosine {
                Some(x) => *x,
                None => -0.001f64,
            };
            let min_magnitude = match min_magnitude {
                Some(x) => *x,
                None => -0.001f64,
            };
            let min_dot_prod = match min_dot_prod {
                Some(x) => *x,
                None => 0i32,
            };

            let records = GenotypeRecords::from_parquet_file(rec);
            assert!(records.is_sorted_by_genome());

            let ind_file = rec.with_extension("ind");
            let _inds = Individuals::from_parquet_file(&ind_file);

            let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);

            // run in parallel and collect row results
            let res: Vec<(u32, u32, i32, i32, i32)> = pairs
                .into_par_iter()
                .map(|(genome1, genome2)| {
                    let mut sum_a2: i32 = 0;
                    let mut sum_b2: i32 = 0;
                    let mut sum_ab: i32 = 0;
                    for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                        let (a, b) = match (a, b) {
                            (Some(a), Some(b)) if a == b => (1, 1),
                            // different rare variants, similarity decrease
                            (Some(_), Some(_)) => (-1, 1),
                            (Some(_), None) => (1, 0),
                            (None, Some(_)) => (0, 1),
                            (None, None) => (0, 0),
                        };
                        sum_a2 += a * a;
                        sum_b2 += b * b;
                        sum_ab += a * b;
                    }

                    let dot_prod = sum_ab;
                    let magnitude = ((sum_a2 as f64) * (sum_b2 as f64)).sqrt();
                    let cosine = dot_prod as f64 / magnitude;

                    let mut out = Some((genome1, genome2, sum_a2, sum_b2, sum_ab));
                    if (dot_prod < min_dot_prod)
                        || (magnitude < min_magnitude)
                        || (cosine < min_cosine)
                    {
                        out = None;
                    }

                    out
                })
                .flatten()
                .collect();

            // for identicial sets, elements correponsing to lower matrix is not
            // calculated use this as an indicator to fill the the low part of
            // the matrix when updating the full jaccard matrix
            let identifical = row_genomes == col_genomes;

            let mut resmat = NamedMatrix::new(row_genomes, col_genomes);

            for (g1, g2, sum_a2, sum_b2, sum_ab) in res {
                let dot_prod = sum_ab;
                let magnitude = ((sum_a2 as f64) * (sum_b2 as f64)).sqrt();
                let cosine = dot_prod as f64 / magnitude;

                if output.is_none() {
                    println!(
                        "g1={g1}, g2={g2}, sum_a2={sum_a2}, sum_b2={sum_b2}, sum_ab={sum_ab}, cosine={cosine:.6}",
                    );
                }

                // update the matrix
                resmat.set_by_names(g1, g2, cosine);
                if identifical {
                    resmat.set_by_names(g2, g1, cosine);
                }
            }

            // write matrix to files
            match output.as_ref() {
                Some(output) => {
                    let dir = output.parent().unwrap().to_str().unwrap();
                    let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                    if !filename.ends_with(".cos") {
                        filename.push_str(".cos");
                    }
                    let mut p = PathBuf::from(dir);
                    p.push(filename);

                    // let resmat0 = resmat.clone();
                    println!("WARN: output option is specified, results are not printed on the screen, check file {:?}", p);
                    // println!("\n writing...");
                    resmat.into_parquet(&p)
                }
                None => {}
            }
        }
        Some(Commands::Grm {
            rec,
            genomes,
            lists,
            min_grm_related,
            output,
        }) => {
            let min_grm_related = match min_grm_related {
                Some(x) => *x,
                None => -0.001f64,
            };

            let mut records = GenotypeRecords::from_parquet_file(rec);
            let ind_file = rec.with_extension("ind");
            let _inds = Individuals::from_parquet_file(&ind_file);
            let sites_file = rec.with_extension("sit");
            let _sites = Sites::from_parquet_file(&sites_file);
            let freq_map = calc_allele_frequency(&mut records, _inds.v().len() * 2, _sites.len());

            records.sort_by_genome();
            assert!(records.is_sorted_by_genome());

            let (pairs, row_genomes, col_genomes) = prep_pairs(&records, genomes, lists);
            let base_sum = calc_base_relationship(&freq_map);

            // run in parallel and collect row results
            let res: Vec<(u32, u32, f64)> = pairs
                .into_par_iter()
                .map(|(genome1, genome2)| {
                    let mut sum = base_sum;
                    for (_pos, a, b) in records.iter_genome_pair_genotypes(genome1, genome2) {
                        let (a, b) = match (a, b) {
                            (Some(a), Some(b)) if a == b => (0.0, 0.0),
                            // different rare variants, similarity decrease
                            (Some(_), Some(_)) => continue, // FIXME
                            (Some(_), None) => (0.0, 2.0),
                            (None, Some(_)) => (2.0, 0.0),
                            (None, None) => continue,
                        };
                        match freq_map.get(&_pos) {
                            Some(p) => {
                                sum -= (2.0 - 2.0 * p) * (2.0 - 2.0 * p) / 2.0 / p / (1.0 - p);
                                sum += (a - 2.0 * p) * (b - 2.0 * p) / 2.0 / p / (1.0 - p);
                            }
                            None => {}
                        }
                    }

                    let n = freq_map.len() as f64;
                    let relationship = sum / n;

                    let mut out = Some((genome1, genome2, relationship));
                    if relationship < min_grm_related {
                        out = None
                    }

                    out
                })
                .flatten()
                .collect();

            // for identicial sets, elements correponsing to lower matrix is not
            // calculated use this as an indicator to fill the the low part of
            // the matrix when updating the full jaccard matrix
            let identifical = row_genomes == col_genomes;

            let mut resmat = NamedMatrix::new(row_genomes, col_genomes);

            for (g1, g2, relationship) in res {
                if output.is_none() {
                    println!("g1={g1}, g2={g2}, relationship={relationship:.6}",);
                }

                // update the matrix
                resmat.set_by_names(g1, g2, relationship);
                if identifical {
                    resmat.set_by_names(g2, g1, relationship);
                }
            }

            // write matrix to files
            match output.as_ref() {
                Some(output) => {
                    let dir = output.parent().unwrap().to_str().unwrap();
                    let mut filename = output.file_name().unwrap().to_str().unwrap().to_owned();
                    if !filename.ends_with(".grm") {
                        filename.push_str(".grm");
                    }
                    let mut p = PathBuf::from(dir);
                    p.push(filename);

                    // let resmat0 = resmat.clone();
                    println!("WARN: output option is specified, results are not printed on the screen, check file {:?}", p);
                    // println!("\n writing...");
                    resmat.into_parquet(&p)
                }
                None => {}
            }
        }
        None => {
            println!("\nUse '-h  or [subcommand] -h' to show help message");
        }
    }
}

// helper functions

fn file_to_u32_vec(p: impl AsRef<Path>) -> Vec<u32> {
    let mut s = String::new();
    let mut reader = std::fs::File::open(p).map(std::io::BufReader::new).unwrap();
    s.clear();
    use std::io::Read;
    reader.read_to_string(&mut s).unwrap();
    let mut v = s
        .trim()
        .split("\n")
        .map(|s| s.parse::<u32>().unwrap())
        .collect::<Vec<u32>>();
    assert!(v.len() != 0);
    v.sort();
    v
}

fn prep_pairs(
    records: &GenotypeRecords,
    genomes: &Option<Vec<u32>>,
    lists: &Vec<PathBuf>,
) -> (
    // pairs
    Vec<(u32, u32)>,
    // row genomes
    Vec<u32>,
    // col genomes
    Vec<u32>,
) {
    // genomes vec
    let genomes: Vec<u32> = match genomes {
        Some(v) => v.clone(),
        None => Vec::new(),
    };

    // read genome lists from files

    let mut lst1 = Vec::<u32>::new();
    let mut lst2 = Vec::<u32>::new();

    match lists.len() {
        0 => {}
        1 => {
            // within subset sharing
            lst1.extend(file_to_u32_vec(&lists[0]));
            lst2.extend_from_slice(lst1.as_slice());
        }
        2 => {
            // inter-subset sharing
            lst1.extend(file_to_u32_vec(&lists[0]));
            lst2.extend(file_to_u32_vec(&lists[1]));

            // ensure non-overlapping
            let mut set = HashSet::<u32>::new();
            set.extend(&lst1);
            let overlapping = lst2.iter().any(|x| set.contains(x));
            if overlapping {
                eprintln!("list1 and list2 are overlapping");
                std::process::exit(-1);
            }
        }
        _ => {
            eprintln!("--lists options can specied no more than 2 times");
            std::process::exit(-1);
        }
    }

    let min_gid = records.records().first().unwrap().get_genome();
    let max_gid = records.records().last().unwrap().get_genome();
    let mut row_genomes = Vec::<u32>::new();
    let mut col_genomes = Vec::<u32>::new();
    let pairs: Vec<(u32, u32)> = match (genomes.len() > 0, lst1.len() > 0) {
        (true, false) => {
            row_genomes.extend(genomes.iter());
            col_genomes.extend(genomes.iter());
            genomes
                .iter()
                .map(|x| *x)
                .cartesian_product(genomes.iter().map(|x| *x))
                .filter(|(a, b)| *a > *b)
                .collect()
        }
        (false, true) => {
            row_genomes.extend(&lst1);
            col_genomes.extend(&lst2);

            if lists.len() == 2 {
                // two non-overallping list. no need to filter
                lst1.iter()
                    .map(|x| *x)
                    .cartesian_product(lst2.iter().map(|x| *x))
                    .collect()
            } else {
                // two identical list need to filtering
                lst1.iter()
                    .map(|x| *x)
                    .cartesian_product(lst2.iter().map(|x| *x))
                    .filter(|(a, b)| *a > *b)
                    .collect()
            }
        }
        (true, true) => {
            panic!("--gnomes and --lists should not be specified at the same time");
        }
        (false, false) => {
            let n = (max_gid - min_gid) as usize;
            if n * n * 8 / 1024 / 1024 / 1024 > 30 {
                eprintln!("too many pairs! consider using --lists option to restrict number of pairs for each run!");
                std::process::exit(-1);
            }
            row_genomes.extend(min_gid..max_gid);
            col_genomes.extend(min_gid..max_gid);
            (min_gid..max_gid)
                .cartesian_product(min_gid..max_gid)
                .filter(|(a, b)| *a > *b)
                .collect()
        }
    };

    (pairs, row_genomes, col_genomes)
}

fn calc_allele_frequency(
    rec: &mut GenotypeRecords,
    num_hap: usize,
    num_sites: usize,
) -> AHashMap<u32, f64> {
    rec.sort_by_position();
    let mut freq_map = AHashMap::<u32, f64>::with_capacity(num_sites);
    let mut pos = u32::MAX;
    let mut cnt = 0u32;
    for r in rec.records() {
        let p = r.get_position();
        if p != pos {
            if pos != u32::MAX {
                freq_map.insert(pos, 1.0 - (cnt as f64) / (num_hap as f64));
            }
            pos = p;
            cnt = 1;
        } else {
            cnt += 1;
        }
        if pos != u32::MAX {
            freq_map.insert(pos, 1.0 - (cnt as f64) / (num_hap as f64));
        }
    }
    freq_map
}

/// Calculate the base sum of relationship (assuming all genotypes are reference allele)
///
/// This is helpful as majority all the genotype are 0s (x=2 references).
/// We can iterative update the run sum of relatedships when we see a non-ref allelel:
/// subtracting a value assume ref allele and add a value consider non-ref
fn calc_base_relationship(freq_map: &AHashMap<u32, f64>) -> f64 {
    freq_map
        .values()
        .map(|p| (2.0 - 2.0 * p) * (2.0 - 2.0 * p) / 2.0 / p / (1.0 - p))
        .sum()
}
