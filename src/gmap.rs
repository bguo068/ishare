pub struct GeneticMap(Vec<(u32, f32)>);
use crate::genome::GenomeInfo;
use csv;
use std::path::Path;

impl GeneticMap {
    pub fn from_genome_info(ginfo: &GenomeInfo) -> Self {
        let mut gw_chr_start_bp = 0u32;
        let mut gw_chr_start_cm = 0.0f32;
        let mut v = Vec::new();

        for (chrlen, plinkmap_fn) in ginfo.chromsize.iter().zip(ginfo.gmaps.iter()) {
            let mut chrmap = GeneticMap::from_plink_map(plinkmap_fn, *chrlen);
            let chrlen_cm = chrmap.get_size_cm();

            chrmap.update_to_genome_wide_coords(gw_chr_start_bp, gw_chr_start_cm);

            v.extend(chrmap.0);

            gw_chr_start_bp += chrlen;
            gw_chr_start_cm += chrlen_cm;
        }
        Self(v)
    }

    pub fn get_gw_chr_start_cm_from_chrid(&self, chrid: usize, ginfo: &GenomeInfo) -> f32 {
        let gw_ch_start_bp = ginfo.gwstarts[chrid];
        self.get_cm(gw_ch_start_bp)
    }

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

    pub fn from_plink_map(p: impl AsRef<Path>, chrlen: u32) -> Self {
        let mut v = vec![(0, 0.0)];
        let mut record = csv::StringRecord::new();

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .from_path(&p)
            .unwrap();

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

        let (x1, y1) = v[v.len() - 2];
        let (x2, y2) = v[v.len() - 1];
        if x2 != chrlen {
            let slope = (y2 - y1) / ((x2 - x1) as f32);
            // println!("slope={slope}");
            let mut end_cm = ((chrlen - x2) as f32) * slope + y2;
            if end_cm < y2 {
                end_cm = y2;
            }
            v.push((chrlen, end_cm));
        }
        // println!("the third last pair: {:?}", v[v.len() - 3]);
        // println!("the second last pair: {:?}", v[v.len() - 2]);
        // println!("last pair: {:?}", v.last().unwrap());
        Self(v)
    }
    pub fn get_cm(&self, bp: u32) -> f32 {
        let idx = self.0.partition_point(|e| e.0 <= bp) - 1;
        let (x1, y1) = self.0[idx];
        let (x2, y2) = self.0[idx + 1];
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
        // println!("0: {:?}", self.0[0]);
        // println!("1: {:?}", self.0[1]);
        // println!("2: {:?}", self.0[2]);
        // println!("idx={idx}, x1={x1}, y1={y1}, x2={x2}, y2={y2}, cm={cm}, slope={slope}, bp={bp}");

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
}
