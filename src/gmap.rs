use crate::genome::GenomeInfo;
use csv;
use serde::{Deserialize, Serialize};
use std::{
    io::{BufWriter, Write},
    path::Path,
};

/// Genetic Map represented as vector of 2-tuple: 0-based bp position and the
/// corresponding cM coordinatesg
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneticMap(Vec<(u32, f32)>);

impl GeneticMap {
    pub fn from_iter(it: impl Iterator<Item = (u32, f32)>) -> Self {
        let v = it.collect();
        Self(v)
    }
    pub fn from_vec(v: Vec<(u32, f32)>) -> Self {
        Self(v)
    }

    pub fn from_constant_rate(rate: f32, chrlen: u32) -> Self {
        Self(vec![(0, 0.0), (chrlen - 1, (chrlen - 1) as f32 / rate)])
    }

    pub fn from_genome_info(ginfo: &GenomeInfo) -> Self {
        let mut gw_chr_start_bp = 0u32;
        let mut gw_chr_start_cm = 0.0f32;
        let mut v = Vec::new();

        for (chrlen, plinkmap_fn) in ginfo.chromsize.iter().zip(ginfo.gmaps.iter()) {
            let mut chrmap = GeneticMap::from_plink_map(plinkmap_fn, *chrlen);
            let chrlen_cm = chrmap.get_size_cm();

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += *chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Self(v)
    }

    pub fn as_slice(&self) -> &[(u32, f32)] {
        &self.0[..]
    }

    /// This method can be used to merge gmaps of all chromosomal gmap into
    /// one for the whole genome
    pub fn from_gmap_vec(gmap_vec: &Vec<Vec<(u32, f32)>>, chromsizes: &Vec<u32>) -> Self {
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
            let (last_bp, last_cm) = v_o_chr.last().unwrap();
            if *last_bp != chrlen - 1 {
                let avg_rate = last_cm / *last_bp as f32;
                v_o_chr.push((chrlen - 1, avg_rate * (chrlen - 1) as f32));
            }

            let mut chrmap = GeneticMap(v_o_chr.clone());
            let chrlen_cm = chrmap.get_size_cm();

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += *chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Self(v)
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
    pub fn from_plink_map(p: impl AsRef<Path>, chrlen: u32) -> Self {
        let mut v = vec![(0, 0.0)];
        let mut record = csv::StringRecord::new();

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(&p)
            .expect(&format!(
                "can not open file {}",
                p.as_ref().to_str().unwrap()
            ));

        while reader.read_record(&mut record).unwrap() {
            // println!("{:?}", record);
            let cm = record[2].parse::<f32>().unwrap();
            // use 0-based position
            let bp: u32 = record[3].parse::<u32>().unwrap() - 1;
            // println!("record: cm={}, bp={}", cm, bp);
            if bp == 0 {
                continue;
            }
            assert!(
                (bp, cm) > *v.last().unwrap(),
                "genetic map should be ordered by position"
            );
            v.push((bp, cm));
        }
        v.sort_by_key(|x| x.0);

        // use average rate for interpolation for the right end; right end
        // should be chrlen-1 not chrlen. Otherwise, there might be collision of
        // bp between end bp of current chromosome with the starting bp of the
        // next chromosome
        let (bp, cm) = v[v.len() - 1];
        let avg_rate = cm / bp as f32;
        assert!(bp <= chrlen, "illegal chromosome end bp in the genetic map");

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
        Self(v)
    }

    pub fn get_cm(&self, bp: u32) -> f32 {
        let idx = self.0.partition_point(|e| e.0 <= bp) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = match self.0.get(idx + 1) {
            Some(x) => *x,
            None => panic!("idx={}, bp={}", idx, bp),
        };
        let slope = (y2 - y1) / (x2 - x1) as f32;
        let mut cm = (bp - x1) as f32 * slope + y1;
        if cm < y1 {
            cm = y1;
        } else if cm > y2 {
            cm = y2;
        } else {
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

    pub fn get_size_cm(&self) -> f32 {
        self.0.last().unwrap().1 - self.0.first().unwrap().1
    }

    pub fn update_to_genome_wide_coords(&mut self, gw_chr_start_bp: u32, gw_chr_start_cm: f32) {
        self.0.iter_mut().for_each(|(bp, cm)| {
            *bp += gw_chr_start_bp;
            *cm += gw_chr_start_cm;
        });
    }

    pub fn to_plink_map_files(&self, ginfo: &GenomeInfo, prefix: impl AsRef<Path>) {
        // make folder
        let parent = prefix.as_ref().parent().unwrap();
        if !parent.exists() {
            std::fs::create_dir_all(parent)
                .expect(&format!("can't creat folder {}", parent.to_string_lossy()));
        }
        let filename = prefix.as_ref().file_name().unwrap().to_str().unwrap();
        for (i, chrname) in ginfo.chromnames.iter().enumerate() {
            // make a file name per chromosome
            let mut path = parent.to_path_buf();
            path.push(format!("{}_{}.map", filename, chrname));
            let mut f = std::fs::File::create(&path).map(BufWriter::new).unwrap();
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

                write!(f, "{} . {} {}\n", chrname, cm, pos).unwrap();
            }
        }
    }
}
