use crate::container::intervals::Intervals;
use crate::gmap::{self, GeneticMap};
use ahash::{HashMap, HashMapExt};
use clap::ValueEnum;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::backtrace::Backtrace;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use snafu::{OptionExt, ResultExt, Snafu};
#[derive(Debug, Snafu)]
pub enum Error {
    IoError {
        source: std::io::Error,
        backtrace: Option<Backtrace>,
    },
    TomlSerError {
        source: Box<toml::ser::Error>,
        backtrace: Option<Backtrace>,
    },
    TomlDeError {
        source: Box<toml::de::Error>,
        backtrace: Option<Backtrace>,
    },
    EmptySliceError {
        backtrace: Option<Backtrace>,
    },
    BincodeError {
        source: bincode::Error,
        backtrace: Option<Backtrace>,
    },
    InferMapFilePrefixError {
        backtrace: Option<Backtrace>,
    },
    #[snafu(transparent)]
    GmapError {
        // use Box here to avoid the circular dependency between two modules
        #[snafu(source(from(gmap::Error, Box::new)))]
        source: Box<gmap::Error>,
    },
}

type Result<T> = std::result::Result<T, Error>;

#[derive(Deserialize, Serialize, Debug)]
struct GenomeFile {
    name: String,
    chromsize: Vec<u32>,
    idx: HashMap<String, usize>,
    chromnames: Vec<String>,
    map_root: Option<PathBuf>,
    gmaps: Vec<PathBuf>,
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

#[derive(Deserialize, Debug, Clone, Serialize, PartialEq, Eq, Default)]
pub struct GenomeInfo {
    pub name: String,
    pub chromsize: Vec<u32>,
    pub chromnames: Vec<String>,
    pub idx: HashMap<String, usize>,
    pub gwstarts: Vec<u32>,
    pub map_root: Option<PathBuf>,
    pub gmaps: Vec<PathBuf>,
}

impl GenomeInfo {
    pub fn new() -> Self {
        Self::default()
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
            map_root: None,
            gmaps: vec![],
        }
    }

    pub fn to_toml_file(&self, p: impl AsRef<Path>) -> Result<()> {
        let gf = GenomeFile {
            name: self.name.clone(),
            chromnames: self.chromnames.clone(),
            chromsize: self.chromsize.clone(),
            map_root: p.as_ref().parent().map(|p| p.to_path_buf()),
            gmaps: self.gmaps.clone(),
            idx: self.idx.clone(),
        };
        let s = toml::to_string(&gf)
            .map_err(Box::new)
            .context(TomlSerSnafu {})?;

        let mut f = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .expect("fail to create file for writing ginfo file");
        f.write_all(s.as_bytes()).context(IoSnafu {})
    }

    pub fn from_toml_file<P>(path: P) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let mut s = String::new();
        let p: &Path = path.as_ref();
        std::fs::File::open(p)
            .context(IoSnafu {})?
            .read_to_string(&mut s)
            .context(IoSnafu {})?;
        let mut gfile: GenomeFile = toml::from_str(&s)
            .map_err(Box::new)
            .context(TomlDeSnafu {})?;
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
            let last = ginfo.gwstarts.last().context(EmptySliceSnafu {})?;
            ginfo.gwstarts.push(*chrlen + *last);
        }
        ginfo.map_root = path.as_ref().parent().map(|p| p.to_path_buf());
        Ok(ginfo)
    }
    pub fn to_gw_pos(&self, chrid: usize, pos: u32) -> u32 {
        self.gwstarts[chrid] + pos
    }
    pub fn to_chr_pos(&self, gw_pos: u32) -> (usize, &str, u32) {
        let chrid = self.gwstarts.partition_point(|x| *x <= (gw_pos)) - 1;
        let pos = gw_pos - self.gwstarts[chrid];
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
            map_root: self.map_root.clone(),
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
                        let mut e = s + parallel_chunksize_bp;
                        if e > chrlen {
                            e = chrlen;
                        }
                        regions.push(Some((chrid, s as u64, Some(e as u64))));
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
    pub fn new_from_name(builtgenome: BuiltinGenome) -> Result<Self> {
        match builtgenome {
            BuiltinGenome::Pf3d7Const15k => Self::new_from_constant_recombination_rate(
                "pf3d7_const15k",
                PF3D7CHRLENS,
                &((1..=14).map(|i| format!("Pf3D7_{i:02}_v3")).collect_vec()),
                0.01 / 15000.0,
            ),
            BuiltinGenome::Sim14chr100cmConst15k => Self::new_from_constant_recombination_rate(
                "sim14chr100cm_const15k",
                &[1_500_000_u32; 14],
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
    ) -> Result<Self> {
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

        let gmap_vec: Vec<_> = chromsizes
            .iter()
            .map(|chrlen| {
                vec![
                    (0, 0.0),
                    (*chrlen - 1, (*chrlen - 1) as f32 * 100.0 * const_recom_rate),
                ]
            })
            .collect();
        let gmap = GeneticMap::from_gmap_vec(&gmap_vec, chromsizes)?;

        Ok(Self { ginfo, gmap })
    }

    pub fn new_from_plink_gmaps(
        genome_name: &str,
        chromsizes: &[u32],
        chromnames: &[String],
        plink_files: &[PathBuf],
    ) -> Result<Self> {
        let mut gfile = GenomeFile {
            name: genome_name.to_owned(),
            chromsize: chromsizes.to_owned(),
            idx: chromnames
                .iter()
                .enumerate()
                .map(|(index, chrname)| (chrname.to_owned(), index))
                .collect(),
            chromnames: chromnames.to_vec(),
            map_root: None,
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
            let last = ginfo.gwstarts.last().context(EmptySliceSnafu {})?;
            ginfo.gwstarts.push(*chrlen + *last);
        }
        let gmap = GeneticMap::from_genome_info(&ginfo)?;
        Ok(Self { ginfo, gmap })
    }

    /// load genome (ginfo and gmap) from a single bincode file, for ease of use
    pub fn load_from_bincode_file(p: impl AsRef<Path>) -> Result<Self> {
        let mut reader = std::fs::File::open(p)
            .map(std::io::BufReader::new)
            .context(IoSnafu {})?;
        bincode::deserialize_from(&mut reader).context(BincodeSnafu {})
    }

    /// save to single bincode file, for ease of use
    pub fn save_to_bincode_file(&self, p: impl AsRef<Path>) -> Result<()> {
        if let Some(parent) = Path::new(p.as_ref()).parent() {
            std::fs::create_dir_all(parent).context(IoSnafu {})?;
        }
        let reader = std::fs::File::create(p)
            .map(std::io::BufWriter::new)
            .context(IoSnafu {})?;
        bincode::serialize_into(reader, &self).context(BincodeSnafu {})
    }

    /// save to readable text file,  including the genome info toml file and a list of per chr recombinaton rate in plink map format
    pub fn save_to_text_files(
        &mut self,
        toml_path: impl AsRef<Path>,
        gmap_path: impl AsRef<Path>,
    ) -> Result<()> {
        if let Some(parent) = toml_path.as_ref().parent() {
            std::fs::create_dir_all(parent).context(IoSnafu {})?;
        }
        self.ginfo.to_toml_file(toml_path.as_ref())?;
        self.set_gmap_path(toml_path, gmap_path);
        self.gmap.to_plink_map_files(&self.ginfo)?;
        Ok(())
    }

    pub fn load_from_text_file(ginfo_path: impl AsRef<Path>) -> Result<Self> {
        let ginfo = GenomeInfo::from_toml_file(ginfo_path)?;
        let gmap = GeneticMap::from_genome_info(&ginfo)?;
        Ok(Self { ginfo, gmap })
    }

    pub fn ginfo(&self) -> &GenomeInfo {
        &self.ginfo
    }

    pub fn gmap(&self) -> &GeneticMap {
        &self.gmap
    }

    pub fn set_gmap_path(&mut self, toml_path: impl AsRef<Path>, gmap_path: impl AsRef<Path>) {
        let dir = &toml_path.as_ref().parent().unwrap();
        self.ginfo.map_root = Some(dir.to_path_buf());
        self.ginfo.gmaps = self
            .ginfo
            .chromnames
            .iter()
            .map(|chrname| {
                dir.join(gmap_path.as_ref())
                    .join(PathBuf::from_str(&format!("{chrname}.map")).unwrap())
            })
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
    let mut genome = Genome::new_from_name(BuiltinGenome::Pf3d7Const15k).unwrap();
    genome.set_gmap_path("tmp_genome.toml", "tmp_gmap");

    genome.save_to_bincode_file("tmp.bincode").unwrap();
    let genome2 = Genome::load_from_bincode_file("tmp.bincode").unwrap();
    assert_eq!(&genome, &genome2);

    genome
        .save_to_text_files("tmp_genome.toml", "tmp_gmap")
        .unwrap();
    let genome3 = Genome::load_from_text_file("tmp_genome.toml").unwrap();

    std::fs::remove_file("tmp_genome.toml").unwrap();
    std::fs::remove_file("tmp.bincode").unwrap();
    std::fs::remove_dir_all("tmp_gmap").unwrap();

    assert_eq!(&genome, &genome3);
}
