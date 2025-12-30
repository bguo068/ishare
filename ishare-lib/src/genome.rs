use crate::container::intervals::Intervals;
use crate::gmap::{self, GeneticMap};
use ahash::{HashMap, HashMapExt};
use bincode::{Decode, Encode};
use clap::ValueEnum;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::backtrace::Backtrace;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use snafu::{OptionExt, ResultExt, Snafu};
#[derive(Debug, Snafu)]
pub enum Error {
    IoError {
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    TomlSerError {
        source: Box<toml::ser::Error>,
        backtrace: Box<Option<Backtrace>>,
    },
    TomlDeError {
        source: Box<toml::de::Error>,
        backtrace: Box<Option<Backtrace>>,
    },
    EmptySliceError {
        backtrace: Box<Option<Backtrace>>,
    },
    BincodeDecodeError {
        #[snafu(source(from(bincode::error::DecodeError, Box::new)))]
        source: Box<bincode::error::DecodeError>,
        backtrace: Box<Option<Backtrace>>,
    },
    BincodeEncodeError {
        #[snafu(source(from(bincode::error::EncodeError, Box::new)))]
        source: Box<bincode::error::EncodeError>,
        backtrace: Box<Option<Backtrace>>,
    },
    InferMapFilePrefixError {
        backtrace: Box<Option<Backtrace>>,
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
    map_root: Option<String>,
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

#[derive(Deserialize, Debug, Encode, Decode, Clone, Serialize, PartialEq, Eq, Default)]
pub struct GenomeInfo {
    pub name: String,
    pub chromsize: Vec<u32>,
    pub chromnames: Vec<String>,
    pub idx: HashMap<String, usize>,
    pub gwstarts: Vec<u32>,
    pub map_root: Option<String>,
    pub gmaps: Vec<String>,
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
            map_root: p
                .as_ref()
                .parent()
                .map(|p| p.to_string_lossy().into_owned()),
            gmaps: self.gmaps.clone(),
            idx: self.idx.clone(),
        };
        let s = toml::to_string(&gf)
            .map_err(Box::new)
            .context(TomlSerSnafu {})?;

        let mut f = std::fs::File::create(p.as_ref())
            .map(BufWriter::new)
            .context(IoSnafu {})?;
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
        ginfo.map_root = path
            .as_ref()
            .parent()
            .map(|p| p.to_string_lossy().into_owned());
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
            let new_chrname = format!("{chrname}_{chrpos}");
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

#[derive(Serialize, Deserialize, Decode, Encode, Debug, PartialEq)]
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
        plink_files: &[String],
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
        bincode::decode_from_reader(&mut reader, bincode::config::standard())
            .context(BincodeDecodeSnafu {})
    }

    /// save to single bincode file, for ease of use
    pub fn save_to_bincode_file(&self, p: impl AsRef<Path>) -> Result<()> {
        if let Some(parent) = Path::new(p.as_ref()).parent() {
            std::fs::create_dir_all(parent).context(IoSnafu {})?;
        }
        let mut writer = std::fs::File::create(p)
            .map(std::io::BufWriter::new)
            .context(IoSnafu {})?;
        bincode::encode_into_std_write(self, &mut writer, bincode::config::standard())
            .context(BincodeEncodeSnafu {})?;
        Ok(())
    }

    /// save to readable text file,  including the genome info toml file and a list of per chr recombinaton rate in plink map format
    pub fn save_to_text_files(&mut self, toml_path: &str, gmap_path: &str) -> Result<()> {
        if let Some(parent) = Path::new(toml_path).parent() {
            std::fs::create_dir_all(parent).context(IoSnafu {})?;
        }
        self.ginfo.to_toml_file(toml_path)?;
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

    pub fn set_gmap_path(&mut self, toml_path: impl AsRef<Path>, gmap_path: &str) {
        let dir = &toml_path.as_ref().parent().unwrap_or(Path::new("."));
        self.ginfo.map_root = Some(dir.to_string_lossy().into_owned());
        self.ginfo.gmaps = self
            .ginfo
            .chromnames
            .iter()
            .map(|chrname| {
                dir.join(gmap_path)
                    .join(format!("{chrname}.map"))
                    .to_string_lossy()
                    .into_owned()
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

#[cfg(test)]
mod comprehensive_tests_group1 {
    use super::*;
    use std::io::Write;

    // TOML parsing robustness tests
    // Temporarily disabled to isolate memory issue
    #[test]
    fn test_malformed_toml_file() {
        let malformed_content = r#"
name = "test_genome"
chromsize = [1000, 2000
# Missing closing bracket and other required fields
"#;

        let temp_file = "test_malformed.toml";
        let mut file = std::fs::File::create(temp_file).unwrap();
        file.write_all(malformed_content.as_bytes()).unwrap();

        let result = GenomeInfo::from_toml_file(temp_file);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), Error::TomlDeError { .. }));

        std::fs::remove_file(temp_file).unwrap();
    }

    #[test]
    fn test_missing_required_fields() {
        let incomplete_content = r#"
name = "test_genome"
# Missing chromsize, chromnames, idx, gmaps
"#;

        let temp_file = "test_incomplete.toml";
        let mut file = std::fs::File::create(temp_file).unwrap();
        file.write_all(incomplete_content.as_bytes()).unwrap();

        let result = GenomeInfo::from_toml_file(temp_file);
        assert!(result.is_err());

        std::fs::remove_file(temp_file).unwrap();
    }

    #[test]
    fn test_invalid_chromosome_sizes() {
        let invalid_content = r#"
name = "test_genome"
chromsize = [1000, 2000]
chromnames = ["chr1", "chr2", "chr3"]
gmaps = ["chr1.map", "chr2.map"]
idx = {"chr1" = 0, "chr2" = 1}
"#;

        let temp_file = "test_invalid_sizes.toml";
        let mut file = std::fs::File::create(temp_file).unwrap();
        file.write_all(invalid_content.as_bytes()).unwrap();

        // This should panic during the check() call due to mismatched lengths
        assert!(std::panic::catch_unwind(|| {
            GenomeInfo::from_toml_file(temp_file).unwrap();
        })
        .is_err());

        std::fs::remove_file(temp_file).unwrap();
    }

    #[test]
    fn test_nonexistent_toml_file() {
        let result = GenomeInfo::from_toml_file("nonexistent_file.toml");
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), Error::IoError { .. }));
    }

    // Coordinate conversion accuracy tests
    #[test]
    fn test_coordinate_conversion_accuracy() {
        let chromsizes = vec![1000, 2000, 1500];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
        let idx: HashMap<String, usize> = [
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
            ("chr3".to_string(), 2),
        ]
        .iter()
        .cloned()
        .collect();
        let gwstarts = vec![0, 1000, 3000];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        // Test genome-wide to chromosome conversion
        let (chrid, chrname, pos) = ginfo.to_chr_pos(0);
        assert_eq!(chrid, 0);
        assert_eq!(chrname, "chr1");
        assert_eq!(pos, 0);

        let (chrid, chrname, pos) = ginfo.to_chr_pos(500);
        assert_eq!(chrid, 0);
        assert_eq!(chrname, "chr1");
        assert_eq!(pos, 500);

        let (chrid, chrname, pos) = ginfo.to_chr_pos(1000);
        assert_eq!(chrid, 1);
        assert_eq!(chrname, "chr2");
        assert_eq!(pos, 0);

        let (chrid, chrname, pos) = ginfo.to_chr_pos(2500);
        assert_eq!(chrid, 1);
        assert_eq!(chrname, "chr2");
        assert_eq!(pos, 1500);

        let (chrid, chrname, pos) = ginfo.to_chr_pos(3000);
        assert_eq!(chrid, 2);
        assert_eq!(chrname, "chr3");
        assert_eq!(pos, 0);

        // Test chromosome to genome-wide conversion
        assert_eq!(ginfo.to_gw_pos(0, 0), 0);
        assert_eq!(ginfo.to_gw_pos(0, 500), 500);
        assert_eq!(ginfo.to_gw_pos(1, 0), 1000);
        assert_eq!(ginfo.to_gw_pos(1, 1500), 2500);
        assert_eq!(ginfo.to_gw_pos(2, 0), 3000);
        assert_eq!(ginfo.to_gw_pos(2, 1000), 4000);
    }

    #[test]
    fn test_chromosome_boundary_conversion() {
        let chromsizes = vec![100, 200];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string()];
        let idx: HashMap<String, usize> = [("chr1".to_string(), 0), ("chr2".to_string(), 1)]
            .iter()
            .cloned()
            .collect();
        let gwstarts = vec![0, 100];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        // Test boundary conditions
        let (chrid, _, pos) = ginfo.to_chr_pos(99);
        assert_eq!(chrid, 0);
        assert_eq!(pos, 99);

        let (chrid, _, pos) = ginfo.to_chr_pos(100);
        assert_eq!(chrid, 1);
        assert_eq!(pos, 0);
    }

    #[test]
    fn test_single_chromosome_conversion() {
        let chromsizes = vec![1000];
        let chromnames = vec!["chr1".to_string()];
        let idx: HashMap<String, usize> = [("chr1".to_string(), 0)].iter().cloned().collect();
        let gwstarts = vec![0];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        let (chrid, chrname, pos) = ginfo.to_chr_pos(500);
        assert_eq!(chrid, 0);
        assert_eq!(chrname, "chr1");
        assert_eq!(pos, 500);

        assert_eq!(ginfo.to_gw_pos(0, 500), 500);
    }

    // Genome partitioning tests
    #[test]
    fn test_partition_genome_with_chunks() {
        let chromsizes = vec![1000, 2000, 500];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
        let idx: HashMap<String, usize> = [
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
            ("chr3".to_string(), 2),
        ]
        .iter()
        .cloned()
        .collect();
        let gwstarts = vec![0, 1000, 3000];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        let partitions = ginfo.partition_genome(Some(600));

        // Should create multiple partitions for larger chromosomes
        assert!(!partitions.is_empty());

        // First chromosome (1000bp) with 600bp chunks should create 2 partitions
        let chr1_partitions: Vec<_> = partitions
            .iter()
            .filter(|p| p.is_some() && p.unwrap().0 == 0)
            .collect();
        assert_eq!(chr1_partitions.len(), 2);

        // Verify partition boundaries
        if let Some(Some((chrid, start, end))) = partitions.first() {
            assert_eq!(*chrid, 0);
            assert_eq!(*start, 0);
            assert!(end.is_some());
        }
    }

    #[test]
    fn test_partition_genome_no_chunking() {
        let chromsizes = vec![1000, 2000];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string()];
        let idx: HashMap<String, usize> = [("chr1".to_string(), 0), ("chr2".to_string(), 1)]
            .iter()
            .cloned()
            .collect();
        let gwstarts = vec![0, 1000];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        let partitions = ginfo.partition_genome(None);

        // Should return single None partition when no chunking
        assert_eq!(partitions.len(), 1);
        assert!(partitions[0].is_none());
    }

    #[test]
    fn test_partition_empty_genome() {
        let ginfo = GenomeInfo::new();
        let partitions = ginfo.partition_genome(Some(1000));
        assert!(partitions.is_empty());
    }

    #[test]
    fn test_partition_small_chunks() {
        let chromsizes = vec![100];
        let chromnames = vec!["chr1".to_string()];
        let idx: HashMap<String, usize> = [("chr1".to_string(), 0)].iter().cloned().collect();
        let gwstarts = vec![0];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        let partitions = ginfo.partition_genome(Some(50));

        // 100bp chromosome with 50bp chunks should create 2 partitions
        assert_eq!(partitions.len(), 2);

        if let Some(Some((chrid, start, end))) = partitions.first() {
            assert_eq!(*chrid, 0);
            assert_eq!(*start, 0);
            assert_eq!(end.unwrap(), 50);
        }

        if let Some(Some((chrid, start, end))) = partitions.get(1) {
            assert_eq!(*chrid, 0);
            assert_eq!(*start, 51);
            assert_eq!(end.unwrap(), 100);
        }
    }
}

#[cfg(test)]
mod comprehensive_tests_group2 {
    use super::*;

    // Built-in genome functionality tests (consolidated to reduce memory usage)
    #[test]
    fn test_builtin_genomes() {
        // Test Pf3d7 genome first
        {
            let genome = Genome::new_from_name(BuiltinGenome::Pf3d7Const15k).unwrap();
            let ginfo = genome.ginfo();

            assert_eq!(ginfo.name, "pf3d7_const15k");
            assert_eq!(ginfo.chromsize.len(), 14);
            assert_eq!(ginfo.chromnames.len(), 14);
            assert_eq!(ginfo.gwstarts.len(), 14);

            // Verify chromosome names format (first few)
            assert_eq!(ginfo.chromnames[0], "Pf3D7_01_v3");
            assert_eq!(ginfo.chromnames[1], "Pf3D7_02_v3");
            assert_eq!(ginfo.chromnames[13], "Pf3D7_14_v3");

            // Verify chromosome sizes match PF3D7CHRLENS (first and last)
            assert_eq!(ginfo.chromsize[0], PF3D7CHRLENS[0]);
            assert_eq!(ginfo.chromsize[13], PF3D7CHRLENS[13]);

            // Verify gwstarts calculation
            assert_eq!(ginfo.gwstarts[0], 0);

            // Verify genetic map exists
            assert!(!genome.gmap().as_slice().is_empty());
        } // Pf3d7 genome dropped here

        // Test simulated genome
        {
            let genome = Genome::new_from_name(BuiltinGenome::Sim14chr100cmConst15k).unwrap();
            let ginfo = genome.ginfo();

            assert_eq!(ginfo.name, "sim14chr100cm_const15k");
            assert_eq!(ginfo.chromsize.len(), 14);
            assert_eq!(ginfo.chromnames.len(), 14);

            // Verify all chromosomes have same size
            for &size in &ginfo.chromsize {
                assert_eq!(size, 1_500_000);
            }

            // Verify chromosome names are numbered 1-14 (check first few)
            assert_eq!(ginfo.chromnames[0], "1");
            assert_eq!(ginfo.chromnames[1], "2");
            assert_eq!(ginfo.chromnames[13], "14");

            // Verify genetic map exists
            assert!(!genome.gmap().as_slice().is_empty());
        } // Sim genome dropped here
    }

    #[test]
    fn test_constant_recombination_rate_validation() {
        let chromsizes = vec![1000, 2000];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string()];
        let const_recom_rate = 0.01 / 15000.0;

        // Create genome in limited scope to ensure cleanup
        {
            let genome = Genome::new_from_constant_recombination_rate(
                "test_genome",
                &chromsizes,
                &chromnames,
                const_recom_rate,
            )
            .unwrap();

            let ginfo = genome.ginfo();
            assert_eq!(ginfo.name, "test_genome");
            assert_eq!(ginfo.chromsize, chromsizes);
            assert_eq!(ginfo.chromnames, chromnames);

            // Verify genetic map has data for 2 chromosomes
            let gmap = genome.gmap();
            assert!(!gmap.as_slice().is_empty());
        } // genome is dropped here
    }

    #[test]
    fn test_total_genome_length() {
        let chromsizes = vec![1000, 2000, 1500];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string(), "chr3".to_string()];
        let idx: HashMap<String, usize> = [
            ("chr1".to_string(), 0),
            ("chr2".to_string(), 1),
            ("chr3".to_string(), 2),
        ]
        .iter()
        .cloned()
        .collect();
        let gwstarts = vec![0, 1000, 3000];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes.clone(),
            chromnames,
            idx,
            gwstarts,
        );

        let total_len = ginfo.get_total_len_bp();
        let expected_len: u32 = chromsizes.iter().sum();
        assert_eq!(total_len, expected_len);
    }
}

#[cfg(test)]
mod comprehensive_tests_group3 {
    use super::*;

    // File I/O operation tests (split to reduce memory usage)
    #[test]
    fn test_bincode_serialization_round_trip() {
        let temp_file = "test_bincode_roundtrip.bincode";

        let genome = Genome::new_from_name(BuiltinGenome::Pf3d7Const15k).unwrap();
        genome.save_to_bincode_file(temp_file).unwrap();
        let loaded_genome = Genome::load_from_bincode_file(temp_file).unwrap();

        assert_eq!(genome, loaded_genome);
        std::fs::remove_file(temp_file).unwrap();
    }
}

#[cfg(test)]
mod comprehensive_tests_group5 {
    use super::*;

    #[test]
    fn test_text_serialization_round_trip() {
        let toml_file = "test_text_roundtrip2.toml";
        let gmap_dir = "test_text_gmap2";

        let mut genome = Genome::new_from_name(BuiltinGenome::Sim14chr100cmConst15k).unwrap();
        genome.set_gmap_path(toml_file, gmap_dir);

        genome.save_to_text_files(toml_file, gmap_dir).unwrap();
        let loaded_genome = Genome::load_from_text_file(toml_file).unwrap();

        assert_eq!(genome.ginfo().name, loaded_genome.ginfo().name);
        assert_eq!(genome.ginfo().chromsize, loaded_genome.ginfo().chromsize);
        assert_eq!(genome.ginfo().chromnames, loaded_genome.ginfo().chromnames);

        std::fs::remove_file(toml_file).unwrap();
        std::fs::remove_dir_all(gmap_dir).unwrap();
    }
}

#[cfg(test)]
mod comprehensive_tests_group6 {
    use super::*;

    #[test]
    fn test_toml_file_creation_with_directories() {
        let deep_path = "test_dir/subdir/genome.toml";
        let ginfo = GenomeInfo::new_from_parts(
            "test".to_string(),
            vec![1000],
            vec!["chr1".to_string()],
            [("chr1".to_string(), 0)].iter().cloned().collect(),
            vec![0],
        );

        // Create parent directories manually
        if let Some(parent) = std::path::Path::new(deep_path).parent() {
            std::fs::create_dir_all(parent).unwrap();
        }

        ginfo.to_toml_file(deep_path).unwrap();

        assert!(std::path::Path::new(deep_path).exists());

        // Cleanup
        std::fs::remove_dir_all("test_dir").unwrap();
    }
}

#[cfg(test)]
mod comprehensive_tests_group7 {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_invalid_bincode_file() {
        let temp_file = "test_invalid.bincode";
        let mut file = std::fs::File::create(temp_file).unwrap();
        // Use a shorter invalid sequence that won't trigger massive allocations
        file.write_all(b"\x00\x01\x02\x03").unwrap();

        let result = Genome::load_from_bincode_file(temp_file);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            Error::BincodeDecodeError { .. }
        ));

        std::fs::remove_file(temp_file).unwrap();
    }

    // Error propagation tests
}

#[cfg(test)]
mod comprehensive_tests_group4 {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_missing_genetic_map_files() {
        let valid_content = r#"
name = "test_genome"
chromsize = [1000, 2000]
chromnames = ["chr1", "chr2"]
gmaps = ["nonexistent1.map", "nonexistent2.map"]
idx = {"chr1" = 0, "chr2" = 1}
"#;

        let temp_file = "test_missing_maps.toml";
        let mut file = std::fs::File::create(temp_file).unwrap();
        file.write_all(valid_content.as_bytes()).unwrap();

        // Loading the genome info should work
        let _ginfo = GenomeInfo::from_toml_file(temp_file).unwrap();

        // But loading a full genome should fail when trying to read genetic maps
        let result = Genome::load_from_text_file(temp_file);
        assert!(result.is_err());

        std::fs::remove_file(temp_file).unwrap();
    }

    #[test]
    fn test_invalid_toml_syntax() {
        let invalid_content = r#"
name = "test_genome"
chromsize = [1000, 2000
chromnames = ["chr1", "chr2"]
# Missing closing bracket
"#;

        let temp_file = "test_invalid_syntax.toml";
        let mut file = std::fs::File::create(temp_file).unwrap();
        file.write_all(invalid_content.as_bytes()).unwrap();

        let result = GenomeInfo::from_toml_file(temp_file);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), Error::TomlDeError { .. }));

        std::fs::remove_file(temp_file).unwrap();
    }

    #[test]
    fn test_genome_from_empty_parts() {
        let ginfo = GenomeInfo::new();

        assert!(ginfo.name.is_empty());
        assert!(ginfo.chromsize.is_empty());
        assert!(ginfo.chromnames.is_empty());
        assert!(ginfo.idx.is_empty());
        assert!(ginfo.gwstarts.is_empty());
        assert!(ginfo.map_root.is_none());
        assert!(ginfo.gmaps.is_empty());
    }

    #[test]
    fn test_split_chromosomes_by_regions() {
        use crate::container::intervals::Intervals;

        // Use smaller data to reduce memory usage
        let chromsizes = vec![100, 200];
        let chromnames = vec!["chr1".to_string(), "chr2".to_string()];
        let idx: HashMap<String, usize> = [("chr1".to_string(), 0), ("chr2".to_string(), 1)]
            .iter()
            .cloned()
            .collect();
        let gwstarts = vec![0, 100];

        let ginfo = GenomeInfo::new_from_parts(
            "test_genome".to_string(),
            chromsizes,
            chromnames,
            idx,
            gwstarts,
        );

        // Create small regions to split on (genome-wide coordinates)
        let regions = Intervals::from_tuples(&[(50, 60), (150, 160)]);

        let split_ginfo = ginfo.split_chromosomes_by_regions(&regions);

        // Should have more chromosomes after splitting
        assert!(split_ginfo.chromnames.len() > ginfo.chromnames.len());
        assert_eq!(split_ginfo.name, "test_genome_rmpeaks");

        // Verify new chromosome names include position info
        assert!(split_ginfo.chromnames.iter().any(|name| name.contains("_")));
    }
}
