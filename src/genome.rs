use crate::container::intervals::Intervals;
use crate::gmap::GeneticMap;
use ahash::{HashMap, HashMapExt};
use clap::ValueEnum;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::io::{BufWriter, Read, Write};
use std::path::Path;

#[derive(Deserialize, Serialize, Debug)]
struct GenomeFile {
    name: String,
    chromsize: Vec<u32>,
    idx: HashMap<String, usize>,
    chromnames: Vec<String>,
    gmaps: Vec<String>,
}

impl GenomeFile {
    fn check(&self) {
        assert_eq!(self.chromsize.len(), self.chromnames.len());
        assert_eq!(self.chromsize.len(), self.gmaps.len());
        assert!(self.chromnames.iter().all(|x| self.idx.contains_key(x)));
    }
}

#[test]
fn test_genome() {
    let mut s = String::new();
    std::fs::File::open("genome.toml")
        .unwrap()
        .read_to_string(&mut s)
        .unwrap();
    let _: GenomeFile = toml::from_str(&s).unwrap();
}

#[derive(Deserialize, Debug, Clone, Serialize, PartialEq, Eq)]
pub struct GenomeInfo {
    pub name: String,
    pub chromsize: Vec<u32>,
    pub chromnames: Vec<String>,
    pub idx: HashMap<String, usize>,
    pub gwstarts: Vec<u32>,
    pub gmaps: Vec<String>,
}

impl GenomeInfo {
    pub fn new() -> Self {
        Self {
            name: String::new(),
            chromsize: vec![],
            chromnames: Vec::new(),
            idx: HashMap::new(),
            gwstarts: vec![],
            gmaps: vec![],
        }
    }

    /// build genomeInfo from exisiting components
    ///
    /// This avoid write to file and then loaded it again to memory
    pub fn new_from_parts(
        name: String,
        chromsize: Vec<u32>,
        chromnames: Vec<String>,
        idx: HashMap<String, usize>,
        gwstarts: Vec<u32>,
    ) -> Self {
        Self {
            name,
            chromsize,
            chromnames,
            idx,
            gwstarts,
            gmaps: vec![],
        }
    }

    pub fn to_toml_file(&self, p: impl AsRef<Path>) {
        let gf = GenomeFile {
            name: self.name.clone(),
            chromnames: self.chromnames.clone(),
            chromsize: self.chromsize.clone(),
            gmaps: self.gmaps.clone(),
            idx: self.idx.clone(),
        };
        let s = toml::to_string(&gf).unwrap();

        let mut f = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .expect("fail to create file for writing ginfo file");
        f.write_all(s.as_bytes()).unwrap();
    }

    pub fn from_toml_file<P>(path: P) -> Self
    where
        P: AsRef<Path>,
    {
        let mut s = String::new();
        let p: &Path = path.as_ref();
        std::fs::File::open(p)
            .unwrap()
            .read_to_string(&mut s)
            .unwrap();
        let mut gfile: GenomeFile = toml::from_str(&s).unwrap();
        gfile.check();
        let mut ginfo = Self::new();
        use std::mem::swap;
        swap(&mut gfile.name, &mut ginfo.name);
        swap(&mut gfile.chromsize, &mut ginfo.chromsize);
        swap(&mut gfile.chromnames, &mut ginfo.chromnames);
        swap(&mut gfile.idx, &mut ginfo.idx);
        swap(&mut gfile.gmaps, &mut ginfo.gmaps);
        ginfo.gwstarts.push(0u32);
        for chrlen in &ginfo.chromsize[0..(ginfo.chromsize.len() - 1)] {
            let last = ginfo.gwstarts.last().unwrap();
            ginfo.gwstarts.push(*chrlen + *last);
        }
        ginfo
    }
    pub fn to_gw_pos(&self, chrid: usize, pos: u32) -> u32 {
        self.gwstarts[chrid] as u32 + pos
    }
    pub fn to_chr_pos(&self, gw_pos: u32) -> (usize, &str, u32) {
        let chrid = self.gwstarts.partition_point(|x| *x <= (gw_pos)) - 1;
        let pos = gw_pos - self.gwstarts[chrid] as u32;
        let chrname = &self.chromnames[chrid];
        (chrid, chrname, pos)
    }
    pub fn get_total_len_bp(&self) -> u32 {
        self.chromsize.iter().sum()
    }

    pub fn split_chromosomes_by_regions(&self, regions: &Intervals<u32>) -> Self {
        let name = format!("{}_rmpeaks", &self.name);
        let mut chromsize = vec![];
        let mut idx = HashMap::new();
        let mut gwstarts = vec![];
        let mut chromnames = vec![];
        let gmaps = vec![];
        let gw_tot_len_bp = self.get_total_len_bp();

        // new gwstart are from old gwstarts plus region boundaries
        gwstarts.extend(self.gwstarts.iter());
        regions
            .iter()
            .for_each(|r| gwstarts.extend(&[r.start, r.end]));
        // dedup
        gwstarts.sort();
        gwstarts.dedup();
        // filter boundaries beyond the total size
        gwstarts.retain(|x| *x < gw_tot_len_bp);

        for (i, gwstart) in gwstarts.iter().enumerate() {
            let (_chrid, chrname, chrpos) = self.to_chr_pos(*gwstart);
            let new_chrname = format!("{}_{}", chrname, chrpos);
            chromnames.push(new_chrname.clone());
            idx.insert(new_chrname, i);
        }

        // use gwstarts to get chrommsome sizes
        gwstarts.push(gw_tot_len_bp);
        chromsize.extend(
            gwstarts
                .iter()
                .zip(gwstarts.iter().skip(1))
                .map(|(a, b)| *b - *a),
        );
        gwstarts.pop();

        Self {
            name,
            chromnames,
            chromsize,
            idx,
            gwstarts,
            gmaps,
        }
    }

    /// return a vector of non-overlapping genome intervals that covers the
    /// whole genomes. This can be use to divide large indexed bcf/vcf files
    /// into smaller chunks for reading in parallel.  
    ///
    /// Each vector element is a 3-tuple, Some(chrid, start_gw_bp,
    /// Option<end_gw_bp>) the end position is can be which means to end of the
    /// correponding chromosome.  if the 3-tuple is a None, it means there is
    /// not segmentation for the the whole genome.
    pub fn partition_genome(
        &self,
        max_len_bp: Option<u32>,
    ) -> Vec<Option<(u32, u64, Option<u64>)>> {
        let mut regions = Vec::new();
        match max_len_bp {
            Some(parallel_chunksize_bp) => {
                for (chrid, chrlen) in self.chromsize.iter().enumerate() {
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
        regions
    }
}

/// A struct that contains both genome info and genetic map
///
/// used to save and read both of genome info and genetic map from from binary or text files

#[derive(Serialize, Deserialize, Debug, PartialEq)]
pub struct Genome {
    ginfo: GenomeInfo,
    gmap: GeneticMap,
}

impl Genome {
    pub fn new_from_name(builtgenome: BuiltinGenome) -> Self {
        match builtgenome {
            BuiltinGenome::Pf3d7Const15k => Self::new_from_constant_recombination_rate(
                "pf3d7_const15k",
                PF3D7CHRLENS,
                &((1..=14).map(|i| format!("Pf3D7_{i:02}_v3")).collect_vec()),
                0.01 / 15000.0,
            ),
            BuiltinGenome::Sim14chr100cmConst15k => Self::new_from_constant_recombination_rate(
                "sim14chr100cm_const15k",
                &[15000_00u32; 14],
                &((1..=14).map(|i| i.to_string()).collect_vec()),
                0.01 / 15000.0,
            ),
        }
    }

    pub fn new_from_constant_recombination_rate(
        genome_name: &str,
        chromsizes: &[u32],
        chromnames: &[String],
        const_recom_rate: f32,
    ) -> Self {
        let m: HashMap<_, _> = chromnames
            .iter()
            .enumerate()
            .map(|(index, chrname)| (chrname.to_owned(), index))
            .collect();

        let mut gw_starts = Vec::with_capacity(chromsizes.len());
        gw_starts.push(0);

        gw_starts.extend(
            chromsizes
                .iter()
                .take(chromsizes.len() - 1)
                .scan(0, |state, &chrlen| {
                    *state += chrlen;
                    Some(*state)
                }),
        );
        let ginfo = GenomeInfo::new_from_parts(
            genome_name.to_owned(),
            chromsizes.to_vec(),
            chromnames.to_vec(),
            m,
            gw_starts.clone(),
        );

        let gmap_vec = chromsizes
            .iter()
            .map(|chrlen| {
                vec![
                    (0, 0.0),
                    (*chrlen - 1, (*chrlen - 1) as f32 * 100.0 * const_recom_rate),
                ]
            })
            .collect();
        let gmap = GeneticMap::from_gmap_vec(&gmap_vec, &chromsizes.to_vec());

        Self { ginfo, gmap }
    }

    pub fn new_from_plink_gmaps(
        genome_name: &str,
        chromsizes: &[u32],
        chromnames: &[String],
        plink_files: &[String],
    ) -> Self {
        let mut gfile = GenomeFile {
            name: genome_name.to_owned(),
            chromsize: chromsizes.to_owned(),
            idx: chromnames
                .iter()
                .enumerate()
                .map(|(index, chrname)| (chrname.to_owned(), index))
                .collect(),
            chromnames: chromnames.to_vec(),
            gmaps: plink_files.to_owned(),
        };
        gfile.check();
        let mut ginfo = GenomeInfo::new();
        use std::mem::swap;
        swap(&mut gfile.name, &mut ginfo.name);
        swap(&mut gfile.chromsize, &mut ginfo.chromsize);
        swap(&mut gfile.chromnames, &mut ginfo.chromnames);
        swap(&mut gfile.idx, &mut ginfo.idx);
        swap(&mut gfile.gmaps, &mut ginfo.gmaps);
        ginfo.gwstarts.push(0u32);
        for chrlen in &ginfo.chromsize[0..(ginfo.chromsize.len() - 1)] {
            let last = ginfo.gwstarts.last().unwrap();
            ginfo.gwstarts.push(*chrlen + *last);
        }
        let gmap = GeneticMap::from_genome_info(&ginfo);
        Self { ginfo, gmap }
    }

    /// load genome (ginfo and gmap) from a single bincode file, for ease of use
    pub fn load_from_bincode_file(p: &str) -> Self {
        let mut reader = std::fs::File::open(p).map(std::io::BufReader::new).unwrap();
        bincode::deserialize_from(&mut reader).unwrap()
    }

    /// save to single bincode file, for ease of use
    pub fn save_to_bincode_file(&self, p: &str) {
        if let Some(parent) = Path::new(p).parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        let reader = std::fs::File::create(p)
            .map(std::io::BufWriter::new)
            .unwrap();
        bincode::serialize_into(reader, &self).unwrap();
    }

    /// save to readable text file,  including the genome info toml file and a list of per chr recombinaton rate in plink map format
    pub fn save_to_text_files(&self, ginfo_path: &str) {
        if let Some(parent) = std::path::Path::new(ginfo_path).parent() {
            std::fs::create_dir_all(parent).unwrap();
        }
        self.ginfo.to_toml_file(ginfo_path);
        let gmap_path_prefix = self.ginfo.gmaps[0]
            .strip_suffix(".map")
            .unwrap()
            .strip_suffix(&self.ginfo.chromnames[0])
            .unwrap()
            .strip_suffix("_")
            .unwrap();
        self.gmap.to_plink_map_files(&self.ginfo, gmap_path_prefix);
    }

    pub fn load_from_text_file(ginfo_path: &str) -> Self {
        let ginfo = GenomeInfo::from_toml_file(ginfo_path);
        let gmap = GeneticMap::from_genome_info(&ginfo);
        Self { ginfo, gmap }
    }

    pub fn ginfo(&self) -> &GenomeInfo {
        &self.ginfo
    }

    pub fn gmap(&self) -> &GeneticMap {
        &self.gmap
    }

    pub fn set_gmap_path_prefix(&mut self, new_gmap_path_prefix: &str) {
        self.ginfo.gmaps = self
            .ginfo
            .chromnames
            .iter()
            .map(|chrname| format!("{new_gmap_path_prefix}_{chrname}.map"))
            .collect::<Vec<_>>();
    }
}

#[derive(Clone, ValueEnum, Debug)]
pub enum BuiltinGenome {
    Pf3d7Const15k,
    Sim14chr100cmConst15k,
}

pub const PF3D7CHRLENS: &[u32] = &[
    640851, 947102, 1067971, 1200490, 1343557, 1418242, 1445207, 1472805, 1541735, 1687656,
    2038340, 2271494, 2925236, 3291936,
];

#[test]
fn test_genome_io() {
    let mut genome = Genome::new_from_name(BuiltinGenome::Pf3d7Const15k);
    genome.set_gmap_path_prefix("tmp_gmap/map");

    genome.save_to_bincode_file("tmp.bincode");
    let genome2 = Genome::load_from_bincode_file("tmp.bincode");
    assert_eq!(&genome, &genome2);

    genome.save_to_text_files("tmp_genome.toml");
    let genome3 = Genome::load_from_text_file("tmp_genome.toml");

    assert_eq!(&genome, &genome3);
}
