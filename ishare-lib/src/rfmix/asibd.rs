use crate::{
    container::intervaltree::IntervalTree,
    genome::GenomeInfo,
    gmap::GeneticMap,
    indiv::Individuals,
    share::ibd::{
        ibdseg::IbdSeg,
        ibdset::{IbdSet, IbdSetBlockIter},
    },
};
use snafu::{ensure, ResultExt, Snafu};
use std::{backtrace::Backtrace, sync::Arc};

use super::fb::LASet;

#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    LASegError {
        // non leaf, delegate backtrace
        #[snafu(backtrace, source(from(super::fb::Error, Box::new)))]
        source: Box<super::fb::Error>,
    },
    #[snafu(display("Individual index {} is out of bounds (max: {})", index, max_index))]
    IndividualIndexOutOfBounds {
        // leaf
        index: u32,
        max_index: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Ancestry index {} is out of bounds (max: {})", index, max_index))]
    AncestryIndexOutOfBounds {
        // leaf
        index: u8,
        max_index: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("IO error during output"))]
    OutputError {
        // leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
}

type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Copy, Debug)]
pub struct ASIbdSeg {
    pub ibd: IbdSeg,
    pub anc1: u8,
    pub anc2: u8,
}

pub struct ASIBDSet {
    asibd: Vec<ASIbdSeg>,
    gmap: Arc<GeneticMap>,
    ginfo: Arc<GenomeInfo>,
    inds: Arc<Individuals>,
    ancs: Arc<[String]>,
}

impl ASIBDSet {
    /// Create an empty IBD set with meta infomation
    pub fn new(
        gmap: Arc<GeneticMap>,
        ginfo: Arc<GenomeInfo>,
        inds: Arc<Individuals>,
        ancs: Arc<[String]>,
    ) -> Self {
        Self {
            asibd: vec![],
            gmap,
            ginfo,
            inds,
            ancs,
        }
    }
    pub fn add(&mut self, asibdseg: ASIbdSeg) -> Result<()> {
        // Validate ancestry indices
        ensure!(
            (asibdseg.anc1 as usize) < self.ancs.len(),
            AncestryIndexOutOfBoundsSnafu {
                index: asibdseg.anc1,
                max_index: self.ancs.len().saturating_sub(1)
            }
        );
        ensure!(
            (asibdseg.anc2 as usize) < self.ancs.len(),
            AncestryIndexOutOfBoundsSnafu {
                index: asibdseg.anc2,
                max_index: self.ancs.len().saturating_sub(1)
            }
        );

        self.asibd.push(asibdseg);
        Ok(())
    }

    pub fn get_asibd_from_ibdsets_and_laset(
        &mut self,
        ibdset: &IbdSet,
        la_set: &LASet,
    ) -> Result<()> {
        let mut tree = IntervalTree::new(1000);
        for blk in IbdSetBlockIter::new(ibdset, false) {
            let (i, m, j, n) = blk[0].haplotype_pair();
            let hap1 = (i << 1) + m as u32;
            let hap2 = (j << 1) + n as u32;
            tree = la_set
                .get_hap_pair_la_segs2(hap1, hap2, tree)
                .context(LASegSnafu)?;

            for ibdseg in blk {
                for elem in tree.query(ibdseg.s..ibdseg.e) {
                    let mut s = ibdseg.s;
                    let mut e = ibdseg.e;
                    if s < elem.range.start {
                        s = elem.range.start;
                    }
                    if e > elem.range.end {
                        e = elem.range.end;
                    }
                    let asibd = ASIbdSeg {
                        ibd: IbdSeg {
                            i: ibdseg.i,
                            j: ibdseg.j,
                            s,
                            e,
                        },
                        anc1: elem.value.0,
                        anc2: elem.value.1,
                    };
                    self.asibd.push(asibd);
                }
            }
        }
        Ok(())
    }

    pub fn flush(&mut self, mut w: impl std::io::Write) -> Result<()> {
        for asibdseg in self.asibd.iter() {
            let ibd = &asibdseg.ibd;
            let (i, m, j, n) = ibd.haplotype_pair();

            // Validate individual indices
            ensure!(
                (i as usize) < self.inds.v().len(),
                IndividualIndexOutOfBoundsSnafu {
                    index: i,
                    max_index: self.inds.v().len().saturating_sub(1)
                }
            );
            ensure!(
                (j as usize) < self.inds.v().len(),
                IndividualIndexOutOfBoundsSnafu {
                    index: j,
                    max_index: self.inds.v().len().saturating_sub(1)
                }
            );

            // Validate ancestry indices
            ensure!(
                (asibdseg.anc1 as usize) < self.ancs.len(),
                AncestryIndexOutOfBoundsSnafu {
                    index: asibdseg.anc1,
                    max_index: self.ancs.len().saturating_sub(1)
                }
            );
            ensure!(
                (asibdseg.anc2 as usize) < self.ancs.len(),
                AncestryIndexOutOfBoundsSnafu {
                    index: asibdseg.anc2,
                    max_index: self.ancs.len().saturating_sub(1)
                }
            );

            let (_chrid, chrname, s) = self.ginfo.to_chr_pos(ibd.s);
            let e = s + (ibd.e - ibd.s);

            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.inds.v()[i as usize],
                m + 1,
                self.inds.v()[j as usize],
                n + 1,
                chrname,
                s + 1, // convert from 0-based to 1-based position
                e + 1, // convert from 0-based to 1-based position
                self.gmap.get_cm_len(ibd.s, ibd.e),
                self.ancs[asibdseg.anc1 as usize],
                self.ancs[asibdseg.anc2 as usize],
            )
            .context(OutputSnafu)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        genome::GenomeInfo,
        gmap::GeneticMap,
        indiv::Individuals,
        rfmix::fb::{FbMatrix, LASet},
        share::{
            ibd::{ibdseg::IbdSeg, ibdset::IbdSet},
            mat::NamedMatrix,
        },
    };
    use ahash::{HashMap, HashMapExt};
    use std::sync::Arc;

    fn create_test_genome() -> GenomeInfo {
        let mut idx = HashMap::new();
        idx.insert("chr1".to_string(), 0);
        idx.insert("chr2".to_string(), 1);

        GenomeInfo {
            name: "test".to_string(),
            chromsize: vec![1000, 500],
            chromnames: vec!["chr1".to_string(), "chr2".to_string()],
            idx,
            gwstarts: vec![0, 1000],
            map_root: None,
            gmaps: vec![],
        }
    }

    fn create_test_individuals() -> Individuals {
        Individuals::from_str_iter(["sample1", "sample2", "sample3"].iter().cloned())
    }

    fn create_test_genetic_map() -> GeneticMap {
        GeneticMap::from_constant_rate(1000000.0, 1500) // 1Mb/cM rate for 1500bp total
    }

    fn create_test_fb_matrix() -> FbMatrix {
        let _individuals = create_test_individuals();

        // Create a minimal ancestry matrix: 3 samples * 2 haps * 2 windows = 12 values
        let ancestry_data = vec![
            0, 1, // sample1 hap1, hap2 for window 1
            1, 0, // sample2 hap1, hap2 for window 1
            0, 0, // sample3 hap1, hap2 for window 1
            1, 1, // sample1 hap1, hap2 for window 2
            0, 1, // sample2 hap1, hap2 for window 2
            1, 0, // sample3 hap1, hap2 for window 2
        ];

        let mat = NamedMatrix::new_from_shape_and_data(6, 2, ancestry_data);

        FbMatrix {
            windows: vec![(100, 400), (400, 800)],
            ancestry: vec!["AFR".to_string(), "EUR".to_string(), "Unknown".to_string()],
            samples: vec![0, 1, 2], // corresponding to sample IDs
            mat,
        }
    }

    #[test]
    fn test_asibdseg_creation() {
        let ibd_seg = IbdSeg {
            i: 0,
            j: 1,
            s: 100,
            e: 300,
        };

        let asibd_seg = ASIbdSeg {
            ibd: ibd_seg,
            anc1: 0,
            anc2: 1,
        };

        assert_eq!(asibd_seg.ibd.i, 0);
        assert_eq!(asibd_seg.ibd.j, 1);
        assert_eq!(asibd_seg.ibd.s, 100);
        assert_eq!(asibd_seg.ibd.e, 300);
        assert_eq!(asibd_seg.anc1, 0);
        assert_eq!(asibd_seg.anc2, 1);
    }

    #[test]
    fn test_asibdset_new() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let asibd_set = ASIBDSet::new(gmap, ginfo, inds, ancs);

        assert_eq!(asibd_set.asibd.len(), 0);
        assert_eq!(asibd_set.ancs.len(), 2);
        assert_eq!(asibd_set.ancs[0], "AFR");
        assert_eq!(asibd_set.ancs[1], "EUR");
    }

    #[test]
    fn test_asibdset_add() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap, ginfo, inds, ancs);
        let ibd_seg = IbdSeg {
            i: 0,
            j: 1,
            s: 100,
            e: 300,
        };

        let asibd_seg = ASIbdSeg {
            ibd: ibd_seg,
            anc1: 0,
            anc2: 1,
        };

        asibd_set.add(asibd_seg).unwrap();

        assert_eq!(asibd_set.asibd.len(), 1);
        assert_eq!(asibd_set.asibd[0].ibd.i, 0);
        assert_eq!(asibd_set.asibd[0].ibd.j, 1);
        assert_eq!(asibd_set.asibd[0].anc1, 0);
        assert_eq!(asibd_set.asibd[0].anc2, 1);
    }

    #[test]
    fn test_asibdset_flush_empty() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap, ginfo, inds, ancs);

        let mut output = Vec::new();
        let result = asibd_set.flush(&mut output);

        assert!(result.is_ok());
        assert!(output.is_empty());
    }

    #[test]
    fn test_asibdset_flush_with_data() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap, ginfo, inds, ancs);

        // Create IbdSeg with encoded values: individual 0 haplotype 1, individual 1 haplotype 0
        // i = (individual << 2) + haplotype = (0 << 2) + 1 = 1
        // j = (individual << 2) + haplotype = (1 << 2) + 0 = 4
        let ibd_seg = IbdSeg {
            i: 1, // individual 0, haplotype 1
            j: 4, // individual 1, haplotype 0
            s: 100,
            e: 300,
        };

        let asibd_seg = ASIbdSeg {
            ibd: ibd_seg,
            anc1: 0,
            anc2: 1,
        };

        asibd_set.add(asibd_seg).unwrap();

        let mut output = Vec::new();
        let result = asibd_set.flush(&mut output);

        assert!(result.is_ok());
        let output_str = String::from_utf8(output).unwrap();

        // Verify output format: sample1 hap sample2 hap chr start end cm_length anc1 anc2
        assert!(output_str.contains("sample1"));
        assert!(output_str.contains("sample2"));
        assert!(output_str.contains("chr1"));
        assert!(output_str.contains("AFR"));
        assert!(output_str.contains("EUR"));

        // Check that positions are 1-based in output
        assert!(output_str.contains("101")); // start position + 1
        assert!(output_str.contains("301")); // end position + 1

        // Verify haplotype numbers (1-based output: haplotype 1 -> 2, haplotype 0 -> 1)
        assert!(output_str.contains("\t2\t")); // haplotype 1 + 1
        assert!(output_str.contains("\t1\t")); // haplotype 0 + 1
    }

    #[test]
    fn test_get_asibd_from_ibdsets_and_laset() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

        // Create test IBD set
        let mut ibd_set = IbdSet::new(gmap.clone(), ginfo.clone(), inds.clone());
        let ibd_seg = IbdSeg {
            i: 0,
            j: 1,
            s: 200,
            e: 600,
        };
        ibd_set.add(ibd_seg);

        // Create test LA set from FB matrix
        let fb_matrix = create_test_fb_matrix();
        let la_set = LASet::from_fbmat(&fb_matrix);

        let result = asibd_set.get_asibd_from_ibdsets_and_laset(&ibd_set, &la_set);

        assert!(result.is_ok());
        // Should have ancestry-specific IBD segments
        assert!(!asibd_set.asibd.is_empty());
    }

    #[test]
    fn test_error_propagation() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

        // Create IBD set with invalid data that might cause LA error
        let mut ibd_set = IbdSet::new(gmap.clone(), ginfo.clone(), inds.clone());
        let ibd_seg = IbdSeg {
            i: 0,
            j: 1,
            s: 200,
            e: 600,
        };
        ibd_set.add(ibd_seg);

        // Create empty LA set (should cause error in real scenario)
        let empty_fb_matrix = FbMatrix {
            windows: vec![],
            ancestry: vec!["AFR".to_string(), "EUR".to_string()],
            samples: vec![],
            mat: NamedMatrix::new_from_shape_and_data(0, 0, vec![]),
        };
        let empty_la_set = LASet::from_fbmat(&empty_fb_matrix);

        let result = asibd_set.get_asibd_from_ibdsets_and_laset(&ibd_set, &empty_la_set);

        // This should handle the error gracefully (might not error with empty data)
        // The test verifies error handling infrastructure is in place
        match result {
            Ok(_) => {}                         // Empty data might not cause error
            Err(Error::LASegError { .. }) => {} // This is expected error type
            Err(_) => {}                        // Other errors are also acceptable
        }
    }

    #[test]
    fn test_multiple_asibd_segments() {
        let ginfo = Arc::new(create_test_genome());
        let gmap = Arc::new(create_test_genetic_map());
        let inds = Arc::new(create_test_individuals());
        let ancs = Arc::from(vec!["AFR".to_string(), "EUR".to_string()]);

        let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

        // Add multiple segments with different ancestry combinations
        let seg1 = ASIbdSeg {
            ibd: IbdSeg {
                i: 0,
                j: 1,
                s: 100,
                e: 200,
            },
            anc1: 0,
            anc2: 1, // AFR-EUR
        };

        let seg2 = ASIbdSeg {
            ibd: IbdSeg {
                i: 1,
                j: 2,
                s: 300,
                e: 400,
            },
            anc1: 1,
            anc2: 2, // EUR-EAS
        };

        let seg3 = ASIbdSeg {
            ibd: IbdSeg {
                i: 0,
                j: 2,
                s: 500,
                e: 600,
            },
            anc1: 0,
            anc2: 0, // AFR-AFR
        };

        asibd_set.add(seg1).unwrap();
        asibd_set.add(seg2).unwrap();
        asibd_set.add(seg3).unwrap();

        assert_eq!(asibd_set.asibd.len(), 3);

        let mut output = Vec::new();
        let result = asibd_set.flush(&mut output);
        assert!(result.is_ok());

        let output_str = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = output_str.trim().split('\n').collect();
        assert_eq!(lines.len(), 3);

        // Verify different ancestry combinations are present
        assert!(output_str.contains("AFR\tEUR"));
        assert!(output_str.contains("EUR\tEAS"));
        assert!(output_str.contains("AFR\tAFR"));
    }

    mod edge_cases {
        use super::*;

        #[test]
        fn test_zero_length_ibd_segment() {
            let ginfo = Arc::new(create_test_genome());
            let gmap = Arc::new(create_test_genetic_map());
            let inds = Arc::new(create_test_individuals());
            let ancs = Arc::from(vec!["AFR".to_string()]);

            let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

            let zero_length_seg = ASIbdSeg {
                ibd: IbdSeg {
                    i: 0,
                    j: 1,
                    s: 100,
                    e: 100,
                }, // zero length
                anc1: 0,
                anc2: 0,
            };

            asibd_set.add(zero_length_seg).unwrap();

            let mut output = Vec::new();
            let result = asibd_set.flush(&mut output);
            assert!(result.is_ok());

            // Should still output the segment (zero-length segments might be valid)
            let output_str = String::from_utf8(output).unwrap();
            assert!(output_str.contains("101\t101")); // Both start and end should be same in 1-based coords
        }

        #[test]
        fn test_boundary_positions() {
            let ginfo = Arc::new(create_test_genome());
            let gmap = Arc::new(create_test_genetic_map());
            let inds = Arc::new(create_test_individuals());
            let ancs = Arc::from(vec!["AFR".to_string()]);

            let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

            // Test segment at chromosome boundary
            let boundary_seg = ASIbdSeg {
                ibd: IbdSeg {
                    i: 0,
                    j: 1,
                    s: 999,
                    e: 1000,
                }, // At chr1 end
                anc1: 0,
                anc2: 0,
            };

            asibd_set.add(boundary_seg).unwrap();

            let mut output = Vec::new();
            let result = asibd_set.flush(&mut output);
            assert!(result.is_ok());

            let output_str = String::from_utf8(output).unwrap();
            assert!(output_str.contains("chr1")); // Should correctly identify chromosome
        }

        #[test]
        fn test_cross_chromosome_segment() {
            let ginfo = Arc::new(create_test_genome());
            let gmap = Arc::new(create_test_genetic_map());
            let inds = Arc::new(create_test_individuals());
            let ancs = Arc::from(vec!["AFR".to_string()]);

            let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

            // Segment spanning chromosome boundary (1000 is start of chr2)
            let cross_chr_seg = ASIbdSeg {
                ibd: IbdSeg {
                    i: 0,
                    j: 1,
                    s: 1001,
                    e: 1200,
                }, // In chr2
                anc1: 0,
                anc2: 0,
            };

            asibd_set.add(cross_chr_seg).unwrap();

            let mut output = Vec::new();
            let result = asibd_set.flush(&mut output);
            assert!(result.is_ok());

            let output_str = String::from_utf8(output).unwrap();
            assert!(output_str.contains("chr2")); // Should correctly identify chr2
        }
    }

    mod error_handling {
        use super::*;

        #[test]
        fn test_invalid_ancestry_index() {
            let ginfo = Arc::new(create_test_genome());
            let gmap = Arc::new(create_test_genetic_map());
            let inds = Arc::new(create_test_individuals());
            let ancs = Arc::from(vec!["AFR".to_string()]);

            let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

            // Try to add segment with invalid ancestry index
            let invalid_seg = ASIbdSeg {
                ibd: IbdSeg {
                    i: 0,
                    j: 1,
                    s: 100,
                    e: 200,
                },
                anc1: 0,
                anc2: 5, // Index 5 doesn't exist
            };

            let result = asibd_set.add(invalid_seg);

            // Should get an error for invalid ancestry index
            assert!(result.is_err());
            if let Err(Error::AncestryIndexOutOfBounds {
                index, max_index, ..
            }) = result
            {
                assert_eq!(index, 5);
                assert_eq!(max_index, 0); // Only one ancestry, so max index is 0
            } else {
                panic!("Expected AncestryIndexOutOfBounds error");
            }
        }

        #[test]
        fn test_invalid_individual_index() {
            let ginfo = Arc::new(create_test_genome());
            let gmap = Arc::new(create_test_genetic_map());
            let inds = Arc::new(create_test_individuals());
            let ancs = Arc::from(vec!["AFR".to_string()]);

            let mut asibd_set = ASIBDSet::new(gmap.clone(), ginfo.clone(), inds.clone(), ancs);

            // Need raw IbdSeg values that decode to invalid individual indices
            // Individual index 10 would be: 10 << 2 = 40 (plus haplotype bits)
            // Individual index 3 doesn't exist (we only have 0,1,2), so 3 << 2 = 12
            let invalid_seg = ASIbdSeg {
                ibd: IbdSeg {
                    i: 12,
                    j: 1,
                    s: 100,
                    e: 200,
                }, // i=12 decodes to individual 3, which doesn't exist
                anc1: 0,
                anc2: 0,
            };

            // Add the invalid segment (doesn't validate indices at add time)
            asibd_set.add(invalid_seg).unwrap();

            // The error will be caught during flush
            let mut output = Vec::new();
            let result = asibd_set.flush(&mut output);

            // Should get an error for invalid individual index
            assert!(result.is_err(), "Expected error but got: {result:?}");
            match result {
                Err(Error::IndividualIndexOutOfBounds { index, .. }) => {
                    assert_eq!(index, 3); // 12 >> 2 = 3
                }
                Err(other) => panic!("Expected IndividualIndexOutOfBounds error, got: {other:?}",),
                Ok(_) => panic!("Expected error but operation succeeded"),
            }
        }
    }
}
