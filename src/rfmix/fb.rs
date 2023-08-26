use itertools::{EitherOrBoth, Itertools};

use crate::{
    container::intervaltree::IntervalTree, genome::GenomeInfo, indiv::Individuals,
    share::mat::NamedMatrix,
};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};
/// rfmix 2 fb.tsv file format
///
/// no. col =  4 + 8 * k * n
///
/// rows 1 columns:
///     - #reference_panel_population:
///     - AFR                
///     - EAS               
///     - EUR                   
///     - NAT
///
/// rows columns:
/// chromosome                    
///        physical_position  
///                  genetic_position  
///                          genetic_marker_index  
///                               8v1_A.NAD_S100:::hap1:::AFR  
///                                        8v1_A.NAD_S100:::hap1:::EAS  
///                                                 8v1_A.NAD_S100:::hap1:::EUR  
///                                                          8v1_A.NAD_S100:::hap1:::NAT  
///                                                                   8v1_A.NAD_S100:::hap2:::AFR  
///                                                                            8v1_A.NAD_S100:::hap2:::EAS  
///                                                                                     8v1_A.NAD_S100:::hap2:::EUR  
///                                                                                              8v1_A.NAD_S100:::hap2:::NAT
/// Following rows:
/// chr22  16747906  2.00077  0   0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16813152  2.16943  5   0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16815536  2.17094  10  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16819340  2.17209  15  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16822150  2.17285  20  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16825214  2.17578  25  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16828406  2.17708  30  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16836097  2.19763  35  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16845542  2.20501  40  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16846476  2.20536  45  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  16934002  2.58794  50  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17044694  2.78064  55  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17131860  3.18325  60  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17132732  3.18399  65  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17134961  3.18621  70  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17139404  3.19298  75  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17143182  3.19878  80  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17148171  3.20484  85  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17151401  3.20783  90  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026
/// chr22  17154637  3.21158  95  0.00000  0.00000  0.00000  1.00000  0.00000  0.00000  0.99974  0.00026

pub struct FbMatrix {
    pub windows: Vec<(u32, u32)>,
    pub ancestry: Vec<String>,
    pub samples: Vec<u32>,
    pub mat: NamedMatrix<u8>,
}

impl FbMatrix {
    pub fn from_fb_csv(
        p: impl AsRef<Path>,
        pos: &[u32],
        ginfo: &GenomeInfo,
        inds: &Individuals,
        min_prob: f32,
        buffer_size_mb: usize,
    ) -> Self {
        let reader = File::open(p.as_ref()).expect(&format!(
            "cannot open file {}",
            p.as_ref().to_str().unwrap()
        ));
        let mut reader = BufReader::with_capacity(1024 * buffer_size_mb, reader);
        let mut buf = String::with_capacity(100000);
        let mut ancestry = Vec::<String>::new();
        let mut windows = vec![];
        let mut win_snp_pos = Vec::<u32>::with_capacity(10000);
        let mut win_snp_idx = Vec::<u32>::with_capacity(10000);
        let mut samples = Vec::<u32>::with_capacity(80);
        let mut mat = Vec::<u8>::with_capacity(1000000);
        let mut v = Vec::<u8>::new();

        // line 1
        let mut ln_cnt = 0;
        reader.read_line(&mut buf).unwrap();
        ln_cnt += 1;
        for a in buf.trim().split("\t").skip(1) {
            ancestry.push(a.to_owned());
        }
        let k_anc = ancestry.len();
        assert!(
            k_anc < u8::MAX as usize,
            "max number of ancestry should be less than 255"
        );
        ancestry.push("Unkown".to_owned());
        buf.clear();

        // line 2, skip 2 columns, read sample name for every k_anc *2 columns
        reader.read_line(&mut buf).unwrap();
        ln_cnt += 1;
        for field in buf.trim().split("\t").skip(4).step_by(k_anc * 2) {
            let sam = field.split(":::").next().unwrap();
            let samid = match inds.m().get(sam) {
                Some(samid) => *samid as u32,
                None => u32::MAX, // if not in individual set to u32::MAX
            };
            samples.push(samid);
        }
        buf.clear();

        // line 3-end
        while reader.read_line(&mut buf).unwrap() > 0 {
            ln_cnt += 1;
            if ln_cnt % 100 == 0 {
                eprint!(
                    "\r\t\tread {:6.3} % positions",
                    100.0 * ln_cnt as f64 / pos.len() as f64
                );
            }
            let mut fields = buf.trim().split("\t");
            let chrname = fields.next().unwrap();
            let chrid = ginfo.idx[chrname];
            let pos: u32 = fields.next().unwrap().parse::<u32>().unwrap() - 1; // parse col 1, use 0-based position
            let gw_pos = ginfo.to_gw_pos(chrid, pos);
            win_snp_pos.push(gw_pos);
            fields.next().unwrap(); // skip col 2
            let _idx: u32 = fields.next().unwrap().parse().unwrap(); // skip col 3
                                                                     // win_snp_idx.push(idx);

            // per site ancestry: Vector
            v.clear();
            v.resize(inds.v().len() * 2, u8::MAX);

            let mut n_chunks = 0;
            let mut n_valid_chunks = 0;
            for chunk in &fields.chunks(k_anc) {
                // this ensures that mat samples order are the same as individuals orders
                let i = n_chunks >> 1;
                let m = n_chunks & 1;
                n_chunks += 1;
                let sid = samples[i];
                if sid == u32::MAX {
                    // skip to next samples if the sample is not in invidual list.
                    continue;
                }
                n_valid_chunks += 1;
                let hapid = (samples[i] << 1) as usize + m;
                //
                let (imax, max) = chunk
                    .enumerate()
                    .map(|(i, p)| (i, p.parse::<f32>().unwrap()))
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .unwrap();
                let anc = match max >= min_prob {
                    true => imax as u8,
                    false => k_anc as u8,
                };
                v[hapid] = anc;
            }
            assert_eq!(n_chunks, samples.len() * 2, "Invalid number of columns!");
            assert_eq!(
                n_valid_chunks,
                inds.v().len() * 2,
                "Invalid number of columns!"
            );
            buf.clear();
            mat.extend(v.iter());
        }

        // calcualte win_snp_idx
        win_snp_idx.extend(
            win_snp_pos
                .into_iter()
                .map(|p| pos.binary_search(&p).unwrap() as u32),
        );

        // caclulate windows
        {
            let mut s = 0;
            for win in 0..win_snp_idx.len() - 1 {
                let e = (win_snp_idx[win] + win_snp_idx[win + 1]) as usize / 2;
                windows.push((pos[s], pos[e]));
                s = e;
            }
            // for last window
            let e = pos.len() - 1;
            windows.push((pos[s], pos[e]));
        }

        let mut mat = NamedMatrix::new_from_shape_and_data(
            win_snp_idx.len() as u32,
            inds.v().len() as u32 * 2,
            mat,
        );
        mat.transpose();

        Self {
            windows,
            samples,
            mat,
            ancestry,
        }
    }

    pub fn get_ancestries(&self) -> &[String] {
        &self.ancestry[..]
    }
}

pub struct LASeg {
    pub win_start: u32,
    pub ancestry: u8,
}

pub struct LASet {
    windows: Vec<(u32, u32)>,
    hap_start_idx: Vec<u32>,
    segs: Vec<LASeg>,
}

impl LASet {
    pub fn from_fbmat(fbmat: &FbMatrix) -> Self {
        let nrows = fbmat.mat.shape().0;
        let mut hap_start_idx = vec![];
        let mut segs = vec![];

        for irow in 0..nrows {
            let hap = fbmat.mat.get_row_slice(irow as u32);
            hap_start_idx.push(segs.len() as u32);

            segs.extend(
                hap.iter()
                    .map(|x| *x)
                    .enumerate()
                    .dedup_by(|a, b| a.1 == b.1)
                    .map(|(w, a)| LASeg {
                        win_start: w as u32,
                        ancestry: a,
                    }),
            );
        }
        Self {
            windows: fbmat.windows.to_owned(),
            hap_start_idx,
            segs,
        }
    }
    pub fn get_lasegs(&self, hap_idx: u32) -> &[LASeg] {
        let hap_idx = hap_idx as usize;
        let s = self.hap_start_idx[hap_idx] as usize;
        let e = match self.hap_start_idx.get(hap_idx + 1) {
            Some(x) => *x as usize,
            None => self.segs.len(),
        };
        &self.segs[s..e]
    }

    /// obtain interval trees for LAsegs for a pair of haplotypes. In each
    /// interval, ancestry is constant for each haplotypes
    ///
    /// Note: buf is needed to calculate end position for each interval
    pub fn get_hap_pair_la_segs(
        &self,
        hap1: u32,
        hap2: u32,
        tree: &mut IntervalTree<u32, (u8, u8)>,
        buf: &mut Vec<(u32, (u8, u8))>,
    ) {
        let segs1 = self.get_lasegs(hap1);
        let segs2 = self.get_lasegs(hap2);
        let mut last_anc1: u8 = 0;
        let mut last_anc2: u8 = 0;
        let iter = segs1
            .iter()
            .merge_join_by(segs2.iter(), |a, b| a.win_start.cmp(&b.win_start))
            .map(|x| match x {
                EitherOrBoth::Both(s1, s2) => {
                    last_anc1 = s1.ancestry;
                    last_anc2 = s2.ancestry;
                    let pos_start = self.windows[s1.win_start as usize].0;
                    (pos_start, (s1.ancestry, s2.ancestry))
                }
                EitherOrBoth::Left(s1) => {
                    last_anc1 = s1.ancestry;
                    let pos_start = self.windows[s1.win_start as usize].0;
                    (pos_start, (s1.ancestry, last_anc2))
                }
                EitherOrBoth::Right(s2) => {
                    last_anc2 = s2.ancestry;
                    let pos_start = self.windows[s2.win_start as usize].0;
                    (pos_start, (last_anc1, s2.ancestry))
                }
            });
        buf.clear();
        buf.extend(iter);
        let last_pos = self.windows.last().unwrap().1;
        buf.push((last_pos, (0, 0)));

        let iter = buf
            .iter()
            .zip(buf.iter().skip(1))
            .map(|(s1, s2)| (s1.0..s2.0, s1.1));

        tree.clear_and_fill_with_iter(iter);
    }

    pub fn get_hap_pair_la_segs2(
        &self,
        hap1: u32,
        hap2: u32,
        mut tree: IntervalTree<u32, (u8, u8)>,
    ) -> IntervalTree<u32, (u8, u8)> {
        let segs1 = self.get_lasegs(hap1);
        let segs2 = self.get_lasegs(hap2);
        let mut last_anc1: u8 = segs1[0].ancestry;
        let mut last_anc2: u8 = segs2[0].ancestry;
        let mut prev_start_pos = self.windows.first().unwrap().0;
        let mut nodes = tree.into_nodes();
        nodes.clear();
        let last_pos = self.windows.last().unwrap().1;
        segs1
            .iter()
            .merge_join_by(segs2.iter(), |a, b| a.win_start.cmp(&b.win_start))
            .skip(1) // skip 1 to avoid add nodes of the first segment twice
            .for_each(|x| match x {
                EitherOrBoth::Both(s1, s2) => {
                    let prev_end_pos = self.windows[s1.win_start as usize].0;
                    let r = prev_start_pos..prev_end_pos;
                    let ancs = (last_anc1, last_anc2);
                    // push previous segment
                    nodes.push(crate::container::intervaltree::Node {
                        element: (r, ancs).into(),
                        max: 0,
                    });
                    prev_start_pos = prev_end_pos;
                    last_anc1 = s1.ancestry;
                    last_anc2 = s2.ancestry;
                }
                EitherOrBoth::Left(s1) => {
                    let prev_end_pos = self.windows[s1.win_start as usize].0;
                    let r = prev_start_pos..prev_end_pos;
                    let ancs = (last_anc1, last_anc2);
                    // push previous segment
                    nodes.push(crate::container::intervaltree::Node {
                        element: (r, ancs).into(),
                        max: 0,
                    });
                    prev_start_pos = prev_end_pos;
                    last_anc1 = s1.ancestry;
                }
                EitherOrBoth::Right(s2) => {
                    let prev_end_pos = self.windows[s2.win_start as usize].0;
                    let r = prev_start_pos..prev_end_pos;
                    let ancs = (last_anc1, last_anc2);
                    // push previous segment
                    nodes.push(crate::container::intervaltree::Node {
                        element: (r, ancs).into(),
                        max: 0,
                    });
                    prev_start_pos = prev_end_pos;
                    last_anc2 = s2.ancestry;
                }
            });
        // push last segment
        nodes.push(crate::container::intervaltree::Node {
            element: ((prev_start_pos..last_pos), (last_anc1, last_anc2)).into(),
            max: 0,
        });

        IntervalTree::new_from_nodes(nodes)
    }
}
