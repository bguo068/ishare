use crate::genome::{self, GenomeInfo};
use bincode::{Decode, Encode};
use csv;
use serde::{Deserialize, Serialize};
use snafu::{ensure, OptionExt, ResultExt, Snafu};
use std::{
    backtrace::Backtrace,
    io::{BufWriter, Write},
    num::{ParseFloatError, ParseIntError},
    path::Path,
};

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    GnomeError {
        #[snafu(source(from(genome::Error, Box::new)))]
        source: Box<genome::Error>,
    },
    #[snafu(display("{source:?}, {path:?}"))]
    CsvError {
        source: csv::Error,
        path: Box<String>,
        backtrace: Box<Option<Backtrace>>,
    },
    IoError {
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    NotEnoughItem {
        backtrace: Box<Option<Backtrace>>,
    },
    MapBpOutOfRange {
        bp: u32,
        chrlen: u32,
        backtrace: Box<Option<Backtrace>>,
    },
    MapItemNotOrderred {
        last_bp: u32,
        bp: u32,
        backtrace: Box<Option<Backtrace>>,
    },
    ParseFloatErr {
        source: ParseFloatError,
        backtrace: Box<Option<Backtrace>>,
    },
    ParseIntErr {
        source: ParseIntError,
        backtrace: Box<Option<Backtrace>>,
    },
    PrefixIsTwoDots {
        backtrace: Box<Option<Backtrace>>,
    },
    PathParentError {
        backtrace: Box<Option<Backtrace>>,
    },
}

type Result<T> = std::result::Result<T, Error>;

/// Genetic Map represented as vector of 2-tuple: 0-based bp position and the
/// corresponding cM coordinatesg
#[derive(Debug, Clone, Serialize, Decode, Encode, Deserialize, PartialEq)]
pub struct GeneticMap(Vec<(u32, f32)>);

impl GeneticMap {
    pub fn from_bp_cm_pair_iter(it: impl Iterator<Item = (u32, f32)>) -> Self {
        let v = it.collect();
        Self(v)
    }
    pub fn from_vec(v: Vec<(u32, f32)>) -> Self {
        Self(v)
    }

    pub fn from_constant_rate(rate: f32, chrlen: u32) -> Self {
        Self(vec![(0, 0.0), (chrlen - 1, (chrlen - 1) as f32 / rate)])
    }

    pub fn from_genome_info(ginfo: &GenomeInfo) -> Result<Self> {
        let mut gw_chr_start_bp = 0u32;
        let mut gw_chr_start_cm = 0.0f32;
        let mut v = Vec::new();

        let dir = ginfo
            .map_root
            .as_ref()
            .map_or(Path::new("."), |p| Path::new(p));
        for (chrlen, plinkmap_fn) in ginfo.chromsize.iter().zip(ginfo.gmaps.iter()) {
            let p = dir.join(plinkmap_fn);
            let mut chrmap = GeneticMap::from_plink_map(p, *chrlen)?;
            let chrlen_cm = chrmap.get_size_cm()?;

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += *chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Ok(Self(v))
    }

    pub fn as_slice(&self) -> &[(u32, f32)] {
        &self.0[..]
    }

    /// This method can be used to merge gmaps of all chromosomal gmap into
    /// one for the whole genome
    pub fn from_gmap_vec(gmap_vec: &[Vec<(u32, f32)>], chromsizes: &[u32]) -> Result<Self> {
        let mut gw_chr_start_bp = 0u32;
        let mut gw_chr_start_cm = 0.0f32;
        let mut v = Vec::new(); // out for a genome
        let mut v_o_chr = Vec::new(); // out for a chromosome

        for (chrlen, gmap_chr) in chromsizes.iter().zip(gmap_vec.iter()) {
            v_o_chr.clear();
            // add left end
            if gmap_chr[0].0 != 0 {
                v_o_chr.push((0, 0.0));
            }
            // add all points
            v_o_chr.extend(gmap_chr.iter());
            // trim on the right size
            v_o_chr.retain(|(bp, _cm)| bp < chrlen);
            // add right end
            let (last_bp, last_cm) = v_o_chr.last().context(NotEnoughItemSnafu {})?;
            if *last_bp != chrlen - 1 {
                let avg_rate = last_cm / *last_bp as f32;
                v_o_chr.push((chrlen - 1, avg_rate * (chrlen - 1) as f32));
            }

            let mut chrmap = GeneticMap(v_o_chr.clone());
            let chrlen_cm = chrmap.get_size_cm()?;

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += *chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Ok(Self(v))
    }

    /// Get chromosome starting cM (in genome-wide space) for a given chromosome
    pub fn get_gw_chr_start_cm_from_chrid(&self, chrid: usize, ginfo: &GenomeInfo) -> f32 {
        let gw_ch_start_bp = ginfo.gwstarts[chrid];
        self.get_cm(gw_ch_start_bp)
    }

    /// Get a vector of chromosome starting cM (in genome-wide space) for all chromosomes
    pub fn get_gw_chr_start_cm_vec(&self, ginfo: &GenomeInfo) -> Vec<f32> {
        let mut v = Vec::new();
        for chrname in ginfo.chromnames.iter() {
            let chrid = ginfo.idx[chrname];
            let gw_ch_start_bp = ginfo.gwstarts[chrid];
            let gw_chr_start_cm = self.get_cm(gw_ch_start_bp);
            v.push(gw_chr_start_cm);
        }

        v
    }

    /// Turn a plink map file into GeneticMap
    ///
    /// Note: currently we require that each plink map should only contain one chromosome
    pub fn from_plink_map(p: impl AsRef<Path>, chrlen: u32) -> Result<Self> {
        let mut v = vec![(0, 0.0)];
        let mut record = csv::StringRecord::new();

        let p_str = p.as_ref().to_string_lossy().into_owned();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(&p)
            .context(CsvSnafu {
                path: p_str.clone(),
            })?;

        while reader.read_record(&mut record).context(CsvSnafu {
            path: p_str.clone(),
        })? {
            // println!("{:?}", record);
            let cm = record[2].parse::<f32>().context(ParseFloatErrSnafu {})?;

            // use 0-based position
            let bp: u32 = record[3].parse::<u32>().context(ParseIntErrSnafu {})? - 1;
            if bp == 0 {
                continue;
            }

            let last_bp = v.last().context(NotEnoughItemSnafu {})?.0;
            ensure!(bp > last_bp, MapItemNotOrderredSnafu { last_bp, bp });
            v.push((bp, cm));
        }
        v.sort_by_key(|x| x.0);

        // use average rate for interpolation for the right end; right end
        // should be chrlen-1 not chrlen. Otherwise, there might be collision of
        // bp between end bp of current chromosome with the starting bp of the
        // next chromosome
        let (bp, cm) = *v.last().context(NotEnoughItemSnafu {})?;
        let avg_rate = cm / bp as f32;
        ensure!(bp <= chrlen, MapBpOutOfRangeSnafu { bp, chrlen });

        let end_bp = chrlen - 1;
        let mut end_cm = ((end_bp - bp) as f32) * avg_rate + cm;
        if end_cm < cm {
            end_cm = cm;
        }
        // fix end bp if the map happend to use chrlen
        if bp == chrlen {
            v.pop();
        }

        // add ends if needed
        if bp != chrlen - 1 {
            v.push((end_bp, end_cm));
        }
        Ok(Self(v))
    }

    pub fn get_cm(&self, bp: u32) -> f32 {
        let idx = self.0.partition_point(|e| e.0 <= bp) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = match self.0.get(idx + 1) {
            Some(x) => *x,
            None => {
                let max_cm = self.get_size_cm().unwrap_or_else(|_| {
                    unreachable!("get_size_cm should not fail as we checked length")
                });
                eprintln!("WARN get_cm: idx={idx}, bp={bp}, cm>={max_cm}");
                return max_cm;
            }
        };
        let slope = (y2 - y1) / (x2 - x1) as f32;
        let mut cm = (bp - x1) as f32 * slope + y1;
        if cm < y1 {
            cm = y1;
        } else if cm > y2 {
            cm = y2;
        }
        cm
    }

    pub fn get_cm_len(&self, s: u32, e: u32) -> f32 {
        self.get_cm(e) - self.get_cm(s)
    }

    pub fn get_bp(&self, cm: f32) -> u32 {
        let idx = self.0.partition_point(|e| e.1 <= cm) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = self.0[idx + 1];
        let slope = (x2 - x1) as f32 / (y2 - y1);
        let mut bp = ((cm - y1) * slope) as u32 + x1;
        if bp < x1 {
            bp = x1;
        } else if bp > x2 {
            bp = x2;
        }
        bp
    }

    pub fn get_size_cm(&self) -> Result<f32> {
        match (self.0.last(), self.0.first()) {
            (Some(last), Some(first)) => Ok(last.1 - first.1),
            _ => NotEnoughItemSnafu.fail(),
        }
    }

    pub fn update_to_genome_wide_coords(&mut self, gw_chr_start_bp: u32, gw_chr_start_cm: f32) {
        self.0.iter_mut().for_each(|(bp, cm)| {
            *bp += gw_chr_start_bp;
            *cm += gw_chr_start_cm;
        });
    }

    pub fn to_plink_map_files(&self, ginfo: &GenomeInfo) -> Result<()> {
        let dir = ginfo
            .map_root
            .as_ref()
            .map_or(Path::new("."), |s| Path::new(s));

        for (i, chrname) in ginfo.chromnames.iter().enumerate() {
            let p = dir.join(&ginfo.gmaps[i]);
            let parent = p.parent().context(PathParentSnafu {})?;
            if !parent.exists() {
                std::fs::create_dir_all(parent).context(IoSnafu {})?;
            }

            let mut f = std::fs::File::create(&p)
                .map(BufWriter::new)
                .context(IoSnafu {})?;
            // find the first end record for each chromosome
            let gwstart = ginfo.gwstarts[i];
            let gwend = match ginfo.gwstarts.get(i + 1) {
                None => ginfo.get_total_len_bp(),
                Some(gwend) => *gwend,
            };
            let s = self.0.partition_point(|x| x.0 < gwstart);
            let e = self.0.partition_point(|x| x.0 < gwend);
            // write records for each chromosomes
            let (pos_offset, cm_offset) = self.0[s];
            for (pos, cm) in &self.0[s..e] {
                let pos = *pos - pos_offset + 1; // 1-based position
                let cm = *cm - cm_offset;

                writeln!(f, "{chrname} . {cm} {pos}").context(IoSnafu {})?;
            }
        }
        Ok(())
    }
}
