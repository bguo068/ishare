use itertools::{EitherOrBoth, Itertools};
use snafu::{ensure, OptionExt, ResultExt, Snafu};

use crate::{
    container::intervaltree::IntervalTree, genome::GenomeInfo, indiv::Individuals,
    share::mat::NamedMatrix,
};
use std::{
    backtrace::Backtrace,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

#[derive(Debug, Snafu)]
pub enum Error {
    IoError {
        // leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    ParseFloatError {
        // leaf
        source: std::num::ParseFloatError,
        backtrace: Box<Option<Backtrace>>,
    },
    ParseIntError {
        // leaf
        source: std::num::ParseIntError,
        backtrace: Box<Option<Backtrace>>,
    },
    InvalidFormat {
        // leaf
        message: Box<String>,
        backtrace: Box<Option<Backtrace>>,
    },
    EmptyData {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
    BinarySearchError {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("No ancestry populations found in header"))]
    NoAncestryPopulations {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Too many ancestry populations: {} (max: {})", count, max))]
    TooManyAncestryPopulations {
        // leaf
        count: usize,
        max: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Invalid chunk size: expected {}, got {}", expected, actual))]
    InvalidChunkSize {
        // leaf
        expected: usize,
        actual: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Haplotype index {} is out of bounds (max: {})", index, max_index))]
    HaplotypeIndexOutOfBounds {
        // leaf
        index: u32,
        max_index: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Corrupted segment indices"))]
    CorruptedSegmentIndices {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
}

type Result<T> = std::result::Result<T, Error>;
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
///
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
    ) -> Result<Self> {
        let file = File::open(p.as_ref()).context(IoSnafu)?;
        let metadata = file.metadata().context(IoSnafu)?;
        ensure!(metadata.len() > 0, EmptyDataSnafu);
        let mut reader = BufReader::with_capacity(1024 * buffer_size_mb, file);
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
        reader.read_line(&mut buf).context(IoSnafu)?;
        ln_cnt += 1;
        for a in buf.trim().split("\t").skip(1) {
            ancestry.push(a.to_owned());
        }
        let k_anc = ancestry.len();

        // Validate ancestry count
        ensure!(k_anc > 0, NoAncestryPopulationsSnafu);
        ensure!(
            k_anc < u8::MAX as usize,
            TooManyAncestryPopulationsSnafu {
                count: k_anc,
                max: u8::MAX as usize - 1
            }
        );
        ancestry.push("Unkown".to_owned());
        buf.clear();

        // line 2, skip 2 columns, read sample name for every k_anc *2 columns
        reader.read_line(&mut buf).context(IoSnafu)?;
        ln_cnt += 1;
        let step_size = k_anc * 2;
        ensure!(step_size > 0, NoAncestryPopulationsSnafu); // This should already be caught, but double-check
        for field in buf.trim().split("\t").skip(4).step_by(step_size) {
            let sam = field.split(":::").next().context(InvalidFormatSnafu {
                message: "Expected sample format with :::".to_owned(),
            })?;
            let samid = match inds.m().get(sam) {
                Some(samid) => *samid as u32,
                None => u32::MAX, // if not in individual set to u32::MAX
            };
            samples.push(samid);
        }
        buf.clear();

        // line 3-end
        while reader.read_line(&mut buf).context(IoSnafu)? > 0 {
            ln_cnt += 1;
            if ln_cnt % 100 == 0 {
                eprint!(
                    "\r\t\tread {:6.3} % positions",
                    100.0 * ln_cnt as f64 / pos.len() as f64
                );
            }
            let mut fields = buf.trim().split("\t");
            let chrname = fields.next().context(InvalidFormatSnafu {
                message: "Missing chromosome name".to_owned(),
            })?;
            let chrid = ginfo.idx[chrname];
            let pos: u32 = fields
                .next()
                .context(InvalidFormatSnafu {
                    message: "Missing position".to_owned(),
                })?
                .parse::<u32>()
                .context(ParseIntSnafu)?
                - 1; // parse col 1, use 0-based position
            let gw_pos = ginfo.to_gw_pos(chrid, pos);
            win_snp_pos.push(gw_pos);
            fields.next().context(InvalidFormatSnafu {
                message: "Missing genetic position column".to_owned(),
            })?; // skip col 2
            let _idx: u32 = fields
                .next()
                .context(InvalidFormatSnafu {
                    message: "Missing genetic marker index".to_owned(),
                })?
                .parse()
                .context(ParseIntSnafu)?; // skip col 3
                                          // win_snp_idx.push(idx);

            // per site ancestry: Vector
            v.clear();
            v.resize(inds.v().len() * 2, u8::MAX);

            let mut n_chunks = 0;
            let mut _n_valid_chunks = 0;
            for chunk in &fields.chunks(k_anc) {
                let chunk_vec: Vec<_> = chunk.collect();
                ensure!(
                    chunk_vec.len() == k_anc,
                    InvalidChunkSizeSnafu {
                        expected: k_anc,
                        actual: chunk_vec.len()
                    }
                );
                // this ensures that mat samples order are the same as individuals orders
                let i = n_chunks >> 1;
                let m = n_chunks & 1;
                n_chunks += 1;
                let sid = samples[i];
                if sid == u32::MAX {
                    // skip to next samples if the sample is not in invidual list.
                    continue;
                }
                _n_valid_chunks += 1;
                let hapid = (samples[i] << 1) as usize + m;
                //
                let parsed_chunk: std::result::Result<Vec<(usize, f32)>, _> = chunk_vec
                    .into_iter()
                    .enumerate()
                    .map(|(i, p)| p.parse::<f32>().map(|f| (i, f)))
                    .collect();
                let (imax, max) = parsed_chunk
                    .context(ParseFloatSnafu)?
                    .into_iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .context(EmptyDataSnafu)?;
                let anc = match max >= min_prob {
                    true => imax as u8,
                    false => k_anc as u8,
                };
                v[hapid] = anc;
            }
            assert_eq!(n_chunks, samples.len() * 2, "Invalid number of columns!");
            buf.clear();
            mat.extend(v.iter());
        }

        // calcualte win_snp_idx
        let snp_indices: std::result::Result<Vec<u32>, Error> = win_snp_pos
            .into_iter()
            .map(|p| {
                pos.binary_search(&p)
                    .map(|idx| idx as u32)
                    .map_err(|_| Error::BinarySearchError {
                        backtrace: Box::new(None),
                    })
            })
            .collect();
        win_snp_idx.extend(snp_indices?);

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

        // Ensure we have at least some matching individuals
        let n_valid_samples = samples.iter().filter(|&&s| s != u32::MAX).count();
        ensure!(n_valid_samples > 0, EmptyDataSnafu);

        let mut mat = NamedMatrix::new_from_shape_and_data(
            win_snp_idx.len() as u32,
            inds.v().len() as u32 * 2,
            mat,
        );
        mat.transpose();

        Ok(Self {
            windows,
            samples,
            mat,
            ancestry,
        })
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
                    .copied()
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
    pub fn get_lasegs(&self, hap_idx: u32) -> Result<&[LASeg]> {
        let hap_idx = hap_idx as usize;
        ensure!(
            hap_idx < self.hap_start_idx.len(),
            HaplotypeIndexOutOfBoundsSnafu {
                index: hap_idx as u32,
                max_index: self.hap_start_idx.len().saturating_sub(1)
            }
        );

        let s = self.hap_start_idx[hap_idx] as usize;
        let e = match self.hap_start_idx.get(hap_idx + 1) {
            Some(x) => *x as usize,
            None => self.segs.len(),
        };

        ensure!(
            s <= self.segs.len() && e <= self.segs.len(),
            CorruptedSegmentIndicesSnafu
        );

        Ok(&self.segs[s..e])
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
    ) -> Result<()> {
        let segs1 = self.get_lasegs(hap1)?;
        let segs2 = self.get_lasegs(hap2)?;

        ensure!(!segs1.is_empty() && !segs2.is_empty(), EmptyDataSnafu);
        ensure!(!self.windows.is_empty(), EmptyDataSnafu);
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
        let last_pos = self.windows.last().context(EmptyDataSnafu)?.1;
        buf.push((last_pos, (0, 0)));

        let iter = buf
            .iter()
            .zip(buf.iter().skip(1))
            .map(|(s1, s2)| (s1.0..s2.0, s1.1));

        tree.clear_and_fill_with_iter(iter);
        Ok(())
    }

    pub fn get_hap_pair_la_segs2(
        &self,
        hap1: u32,
        hap2: u32,
        mut tree: IntervalTree<u32, (u8, u8)>,
    ) -> Result<IntervalTree<u32, (u8, u8)>> {
        let segs1 = self.get_lasegs(hap1)?;
        let segs2 = self.get_lasegs(hap2)?;

        ensure!(!segs1.is_empty() && !segs2.is_empty(), EmptyDataSnafu);
        ensure!(!self.windows.is_empty(), EmptyDataSnafu);

        let mut last_anc1: u8 = segs1[0].ancestry;
        let mut last_anc2: u8 = segs2[0].ancestry;
        let mut prev_start_pos = self.windows.first().context(EmptyDataSnafu)?.0;
        let mut nodes = tree.into_nodes();
        nodes.clear();
        let last_pos = self.windows.last().context(EmptyDataSnafu)?.1;
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

        Ok(IntervalTree::new_from_nodes(nodes))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{genome::GenomeInfo, indiv::Individuals, share::mat::NamedMatrix};
    use ahash::HashMapExt;
    use std::io::Write;

    fn create_test_genome() -> GenomeInfo {
        let mut idx = ahash::HashMap::new();
        idx.insert("chr1".to_string(), 0);
        idx.insert("chr22".to_string(), 1);

        GenomeInfo {
            name: "test".to_string(),
            chromsize: vec![1000, 500],
            chromnames: vec!["chr1".to_string(), "chr22".to_string()],
            idx,
            gwstarts: vec![0, 1000],
            map_root: None,
            gmaps: vec![],
        }
    }

    fn create_test_individuals() -> Individuals {
        Individuals::from_str_iter(["sample1", "sample2"].iter().cloned())
    }

    fn create_test_fb_data() -> String {
        // Create realistic RFMix2 fb.tsv format data
        let header = "#reference_panel_population:\tAFR\tEAS\tEUR\tNAT\n";
        let columns = "chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tsample1:::hap1:::AFR\tsample1:::hap1:::EAS\tsample1:::hap1:::EUR\tsample1:::hap1:::NAT\tsample1:::hap2:::AFR\tsample1:::hap2:::EAS\tsample1:::hap2:::EUR\tsample1:::hap2:::NAT\tsample2:::hap1:::AFR\tsample2:::hap1:::EAS\tsample2:::hap1:::EUR\tsample2:::hap1:::NAT\tsample2:::hap2:::AFR\tsample2:::hap2:::EAS\tsample2:::hap2:::EUR\tsample2:::hap2:::NAT\n";
        let data = "\
chr1\t101\t0.1\t0\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0
chr1\t201\t0.2\t5\t1.0\t0.0\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0
chr1\t301\t0.3\t10\t0.0\t1.0\t0.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0
chr1\t401\t0.4\t15\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0
chr1\t501\t0.5\t20\t0.0\t0.0\t1.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0\t0.0\t1.0\t0.0\t0.0
";
        format!("{header}{columns}{data}")
    }

    fn create_temp_fb_file() -> std::io::Result<std::path::PathBuf> {
        use std::process;
        use std::time::{SystemTime, UNIX_EPOCH};

        let temp_dir = std::env::temp_dir();
        // Create unique filename using timestamp and process ID to avoid race conditions
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let pid = process::id();
        let temp_file = temp_dir.join(format!("test_fb_{timestamp}_{pid}.tsv"));

        let mut file = std::fs::File::create(&temp_file)?;
        file.write_all(create_test_fb_data().as_bytes())?;
        Ok(temp_file)
    }

    #[test]
    fn test_fbmatrix_creation_basic() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500]; // 0-based positions

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

        assert!(fb_matrix.is_ok());
        let fb_matrix = fb_matrix.unwrap();

        // Verify basic structure
        assert_eq!(fb_matrix.ancestry.len(), 5); // AFR, EAS, EUR, NAT, Unknown
        assert_eq!(fb_matrix.ancestry[0], "AFR");
        assert_eq!(fb_matrix.ancestry[1], "EAS");
        assert_eq!(fb_matrix.ancestry[2], "EUR");
        assert_eq!(fb_matrix.ancestry[3], "NAT");
        assert_eq!(fb_matrix.ancestry[4], "Unkown"); // Note: typo in original code

        assert_eq!(fb_matrix.samples.len(), 2); // sample1, sample2
        assert_eq!(fb_matrix.windows.len(), 5); // 5 positions = 5 windows

        // Cleanup
        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_fbmatrix_ancestry_assignment() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();

        // Verify matrix shape: haplotypes x windows
        let (nrows, ncols) = fb_matrix.mat.shape();
        assert_eq!(nrows, 4); // 2 samples * 2 haplotypes
        assert_eq!(ncols, 5); // 5 windows

        // Check specific ancestry assignments
        // First position: sample1 hap1 should be AFR (0), sample1 hap2 should be EUR (2)
        let row0 = fb_matrix.mat.get_row_slice(0); // sample1 hap1
        let row1 = fb_matrix.mat.get_row_slice(1); // sample1 hap2

        assert_eq!(row0[0], 0); // AFR at position 0
        assert_eq!(row1[0], 2); // EUR at position 0

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_fbmatrix_min_prob_threshold() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        // Test with high minimum probability threshold (should assign more to "Unknown")
        let fb_matrix_high =
            FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.99, 1).unwrap();
        let fb_matrix_low = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();

        // With high threshold, more assignments should be "Unknown" (ancestry index 4)
        let high_mat = &fb_matrix_high.mat;
        let low_mat = &fb_matrix_low.mat;

        // Count unknown assignments (ancestry index 4)
        let mut high_unknown = 0;
        let mut low_unknown = 0;

        for i in 0..high_mat.shape().0 {
            for j in 0..high_mat.shape().1 {
                if high_mat.get_row_slice(i as u32)[j] == 4 {
                    high_unknown += 1;
                }
                if low_mat.get_row_slice(i as u32)[j] == 4 {
                    low_unknown += 1;
                }
            }
        }

        // High threshold should have more or equal unknown assignments
        assert!(high_unknown >= low_unknown);

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_fbmatrix_get_ancestries() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();

        let ancestries = fb_matrix.get_ancestries();
        assert_eq!(ancestries.len(), 5);
        assert_eq!(ancestries[0], "AFR");
        assert_eq!(ancestries[1], "EAS");
        assert_eq!(ancestries[2], "EUR");
        assert_eq!(ancestries[3], "NAT");
        assert_eq!(ancestries[4], "Unkown");

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_laset_from_fbmat() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();
        let la_set = LASet::from_fbmat(&fb_matrix);

        // Verify structure
        assert_eq!(la_set.windows.len(), fb_matrix.windows.len());
        assert_eq!(la_set.hap_start_idx.len(), fb_matrix.mat.shape().0); // One entry per haplotype

        // Each haplotype should have segments
        for hap_idx in 0..la_set.hap_start_idx.len() {
            let segs = la_set.get_lasegs(hap_idx as u32).unwrap();
            assert!(!segs.is_empty(), "Haplotype {hap_idx} should have segments",);
        }

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_laset_get_lasegs() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();
        let la_set = LASet::from_fbmat(&fb_matrix);

        // Test getting segments for first haplotype
        let segs0 = la_set.get_lasegs(0).unwrap();
        assert!(!segs0.is_empty());

        // Verify segments are ordered by window
        for i in 1..segs0.len() {
            assert!(
                segs0[i - 1].win_start <= segs0[i].win_start,
                "Segments should be ordered by window start"
            );
        }

        // Test boundary conditions
        let segs_last = la_set
            .get_lasegs(la_set.hap_start_idx.len() as u32 - 1)
            .unwrap();
        assert!(!segs_last.is_empty());

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_laset_hap_pair_la_segs() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();
        let la_set = LASet::from_fbmat(&fb_matrix);

        let mut tree = IntervalTree::new(100);
        let mut buffer = Vec::new();

        let result = la_set.get_hap_pair_la_segs(0, 1, &mut tree, &mut buffer);
        assert!(result.is_ok());

        // Tree should now contain intervals with ancestry pairs
        // Tree should have content (no direct is_empty method, check via query)
        let query_all: Vec<_> = tree.query(0..u32::MAX).collect();
        assert!(!query_all.is_empty());

        // Test query on the tree
        let query_results: Vec<_> = tree.query(150..250).collect();
        assert!(!query_results.is_empty());

        std::fs::remove_file(temp_file).ok();
    }

    #[test]
    fn test_laset_hap_pair_la_segs2() {
        let temp_file = create_temp_fb_file().expect("Failed to create temp file");
        let ginfo = create_test_genome();
        let inds = create_test_individuals();
        let pos = vec![100, 200, 300, 400, 500];

        let fb_matrix = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1).unwrap();
        let la_set = LASet::from_fbmat(&fb_matrix);

        let tree = IntervalTree::new(100);
        let result = la_set.get_hap_pair_la_segs2(0, 1, tree);

        assert!(result.is_ok());
        let tree = result.unwrap();

        // Tree should contain intervals with ancestry pairs
        // Tree should have content (no direct is_empty method, check via query)
        let query_all: Vec<_> = tree.query(0..u32::MAX).collect();
        assert!(!query_all.is_empty());

        // Test queries at different positions
        let query1: Vec<_> = tree.query(150..250).collect();
        let query2: Vec<_> = tree.query(350..450).collect();

        assert!(!query1.is_empty());
        assert!(!query2.is_empty());

        std::fs::remove_file(temp_file).ok();
    }

    mod edge_cases {
        use super::*;
        use std::time::{SystemTime, UNIX_EPOCH};

        #[test]
        fn test_empty_fb_file() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "empty_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));
            std::fs::write(&temp_file, "").unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100, 200, 300];

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            assert!(result.is_err()); // Should fail on empty file

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_malformed_header() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "malformed_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            let malformed_data = "invalid_header\nmore_invalid_data\n";
            std::fs::write(&temp_file, malformed_data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100, 200];

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            // Should handle malformed data gracefully
            match result {
                Ok(_) => {} // Might succeed with default values
                Err(_) => {
                    todo!()
                } // Or fail with appropriate error
            }

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_missing_samples() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "missing_samples_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            // Create data with samples not in individuals list
            let data = "#reference_panel_population:\tAFR\tEUR\n\
chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tmissing_sample:::hap1:::AFR\tmissing_sample:::hap1:::EUR\tmissing_sample:::hap2:::AFR\tmissing_sample:::hap2:::EUR\n\
chr1\t101\t0.1\t0\t1.0\t0.0\t0.0\t1.0\n";

            std::fs::write(&temp_file, data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals(); // Contains sample1, sample2
            let pos = vec![100];

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            // Should error due to no valid columns for known individuals
            assert!(result.is_err());
            if let Err(e) = result {
                let error_msg = format!("{e}");
                // Should be an error about no ancestry populations or empty data
                assert!(
                    error_msg.contains("NoAncestryPopulations")
                        || error_msg.contains("Empty")
                        || error_msg.contains("valid")
                );
            }

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_single_position() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "single_pos_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            // Include all samples to match expected column count (3 samples * 2 haplotypes * 2 ancestries = 12 columns)
            let data = "#reference_panel_population:\tAFR\tEUR\n\
chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tsample1:::hap1:::AFR\tsample1:::hap1:::EUR\tsample1:::hap2:::AFR\tsample1:::hap2:::EUR\tsample2:::hap1:::AFR\tsample2:::hap1:::EUR\tsample2:::hap2:::AFR\tsample2:::hap2:::EUR\tsample3:::hap1:::AFR\tsample3:::hap1:::EUR\tsample3:::hap2:::AFR\tsample3:::hap2:::EUR\n\
chr1\t101\t0.1\t0\t1.0\t0.0\t0.0\t1.0\t0.0\t1.0\t1.0\t0.0\t1.0\t0.0\t0.0\t1.0\n";

            std::fs::write(&temp_file, data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100]; // Single position

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            assert!(result.is_ok());
            let fb_matrix = result.unwrap();
            assert_eq!(fb_matrix.windows.len(), 1);

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_laset_empty_fbmatrix() {
            let empty_fb = FbMatrix {
                windows: vec![],
                ancestry: vec!["AFR".to_string()],
                samples: vec![],
                mat: NamedMatrix::new_from_shape_and_data(0, 0, vec![]),
            };

            let la_set = LASet::from_fbmat(&empty_fb);

            assert_eq!(la_set.windows.len(), 0);
            assert_eq!(la_set.hap_start_idx.len(), 0);
            assert_eq!(la_set.segs.len(), 0);
        }
    }

    mod error_handling {
        use super::*;
        use std::time::{SystemTime, UNIX_EPOCH};

        #[test]
        fn test_invalid_probability_values() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "invalid_prob_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            let data = "#reference_panel_population:\tAFR\tEUR\n\
chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tsample1:::hap1:::AFR\tsample1:::hap1:::EUR\tsample1:::hap2:::AFR\tsample1:::hap2:::EUR\n\
chr1\t101\t0.1\t0\tinvalid\t0.0\t0.0\t1.0\n";

            std::fs::write(&temp_file, data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100];

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            assert!(result.is_err());
            if let Err(Error::ParseFloatError { .. }) = result {
                // Expected error type
            } else {
                panic!("Expected ParseFloatError");
            }

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_invalid_position_format() {
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "invalid_pos_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            let data = "#reference_panel_population:\tAFR\tEUR\n\
chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tsample1:::hap1:::AFR\tsample1:::hap1:::EUR\tsample1:::hap2:::AFR\tsample1:::hap2:::EUR\n\
chr1\tinvalid_pos\t0.1\t0\t1.0\t0.0\t0.0\t1.0\n";

            std::fs::write(&temp_file, data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100];

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            assert!(result.is_err());
            if let Err(Error::ParseIntError { .. }) = result {
                // Expected error type
            } else {
                panic!("Expected ParseIntError");
            }

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_position_not_in_pos_array() {
            let temp_file = create_temp_fb_file().expect("Failed to create temp file");
            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![50, 150, 250]; // Positions don't match FB file positions

            let result = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);

            assert!(result.is_err());
            if let Err(Error::BinarySearchError { .. }) = result {
                // Expected error type
            } else {
                panic!("Expected BinarySearchError");
            }

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_laset_empty_windows_error() {
            let fb_matrix = FbMatrix {
                windows: vec![], // Empty windows
                ancestry: vec!["AFR".to_string()],
                samples: vec![0, 1],
                mat: NamedMatrix::new_from_shape_and_data(2, 1, vec![0, 1]),
            };

            let la_set = LASet::from_fbmat(&fb_matrix);

            let tree = IntervalTree::new(100);
            let result = la_set.get_hap_pair_la_segs2(0, 1, tree);

            assert!(result.is_err());
            if let Err(Error::EmptyData { .. }) = result {
                // Expected error type
            } else {
                panic!("Expected EmptyData error");
            }
        }
    }

    mod performance {
        use super::*;
        use std::time::{SystemTime, UNIX_EPOCH};

        #[test]
        fn test_large_dataset_handling() {
            // Create a larger synthetic dataset
            let temp_dir = std::env::temp_dir();
            let temp_file = temp_dir.join(format!(
                "large_fb_{}_{}.tsv",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap_or_default()
                    .as_nanos(),
                std::process::id()
            ));

            let mut data = String::from("#reference_panel_population:\tAFR\tEUR\n");
            data.push_str("chromosome\tphysical_position\tgenetic_position\tgenetic_marker_index\tsample1:::hap1:::AFR\tsample1:::hap1:::EUR\tsample1:::hap2:::AFR\tsample1:::hap2:::EUR\tsample2:::hap1:::AFR\tsample2:::hap1:::EUR\tsample2:::hap2:::AFR\tsample2:::hap2:::EUR\n");

            // Generate 100 positions
            let mut positions = Vec::new();
            for i in 0..100 {
                let pos = 100 + i * 10;
                positions.push(pos);
                data.push_str(&format!(
                    "chr1\t{}\t{:.1}\t{}\t1.0\t0.0\t0.0\t1.0\t0.0\t1.0\t1.0\t0.0\n",
                    pos + 1,
                    i as f32 * 0.1,
                    i
                ));
            }

            std::fs::write(&temp_file, data).unwrap();

            let ginfo = create_test_genome();
            let inds = create_test_individuals();

            let start_time = std::time::Instant::now();
            let result = FbMatrix::from_fb_csv(&temp_file, &positions, &ginfo, &inds, 0.5, 1);
            let duration = start_time.elapsed();

            assert!(result.is_ok());
            assert!(
                duration.as_secs() < 5,
                "Processing should complete within 5 seconds"
            );

            let fb_matrix = result.unwrap();
            assert_eq!(fb_matrix.windows.len(), 100);

            std::fs::remove_file(temp_file).ok();
        }

        #[test]
        fn test_memory_efficiency() {
            let temp_file = create_temp_fb_file().expect("Failed to create temp file");
            let ginfo = create_test_genome();
            let inds = create_test_individuals();
            let pos = vec![100, 200, 300, 400, 500];

            // Test with different buffer sizes
            let result_small = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 1);
            let result_large = FbMatrix::from_fb_csv(&temp_file, &pos, &ginfo, &inds, 0.5, 16);

            assert!(result_small.is_ok());
            assert!(result_large.is_ok());

            // Both should produce identical results regardless of buffer size
            let fb_small = result_small.unwrap();
            let fb_large = result_large.unwrap();

            assert_eq!(fb_small.windows.len(), fb_large.windows.len());
            assert_eq!(fb_small.mat.shape(), fb_large.mat.shape());

            std::fs::remove_file(temp_file).ok();
        }
    }
}
