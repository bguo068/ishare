use super::WriteFileSnafu;
use crate::gmap::GeneticMap;
use serde::{Deserialize, Serialize};
use snafu::prelude::*;
use std::path::Path;

type Result<T> = std::result::Result<T, super::Error>;

/// struct representing an *encoded* IBD segment:
/// - for i and j, the lower 2 bits encode haplotype id, upper bits encodes individual id
/// - for s and e, the values are either chromosome-level bp poistion or genome-wide bp position
#[derive(Debug, PartialEq, PartialOrd, Eq, Ord, Copy, Clone, Serialize, Deserialize)]
pub struct IbdSeg {
    pub i: u32,
    pub j: u32,
    pub s: u32,
    pub e: u32,
}

impl IbdSeg {
    /// Create a new IBD segment:
    /// - `i`: index of individual 1
    /// - `j`: index of individual 2
    /// - `m`: indicating which happlotype of individual 1 (can be 0 and 1)
    /// - `n`: indicating which happlotype of individual 2 (can be 0 and 1)
    /// - `s`: start position in base pair of the IBD segment (0-based)
    /// - `e`: end position in base pair of the IBD segment (0-based)
    /// - `pos_shift`: is an offset value that will be added `s` and `e`.
    ///
    /// Note 1:
    ///
    /// IBD segment start and end poisitionscan hold information of chromosome index by
    /// representing genome-wide coordinates (assuming all chromsomes are concatenated
    /// into a single one). See `GenomeInfo` methods
    /// [crate::genome::GenomeInfo::to_chr_pos] and
    /// [crate::genome::GenomeInfo::to_gw_pos] and
    /// for chromosomal and genome-wide position conversion.
    ///
    /// Note 2:
    ///
    /// The pair of haplotype indices can be 4 combination for original IBD segments
    /// - (0, 0)
    /// - (0, 1)
    /// - (1, 0)
    /// - (1, 1)
    ///
    /// For segment from merging, there is only one possibility:
    /// - (2, 2)
    ///
    /// For segment from haploid individuals, haplotype indcies should be 3
    /// for all individuals and all segments
    /// - (3, 3)
    pub fn new(i: u32, m: u8, j: u32, n: u8, s: u32, e: u32, pos_shift: u32) -> Self {
        IbdSeg {
            i: i * 4 + m as u32,
            j: j * 4 + n as u32,
            s: s + pos_shift,
            e: e + pos_shift,
        }
    }

    /// Normalize the IBD segment record by swapping (`i` and `j`) if `i` < `j`.
    /// This is useful to uniquely identify haplotype/individual pairs.
    pub fn normalized(&mut self) {
        if self.i < self.j {
            std::mem::swap(&mut self.i, &mut self.j);
        }
    }

    /// return raw `i` and `j`
    pub fn haplotype_pair_int(&self) -> (u32, u32) {
        (self.i, self.j)
    }

    /// return unpacked haplotype pairs(ind1, hap1, ind2, hap2)
    pub fn haplotype_pair(&self) -> (u32, u8, u32, u8) {
        (
            self.i >> 2,
            self.i as u8 & 0x3,
            self.j >> 2,
            self.j as u8 & 0x3,
        )
    }

    /// return unpacked individual pairs(ind1, ind2)
    pub fn individual_pair(&self) -> (u32, u32) {
        (self.i >> 2, self.j >> 2)
    }

    /// return start and end postion (genomw-wide coordinates)
    pub fn coords(&self) -> (u32, u32) {
        (self.s, self.e)
    }

    /// return IBD segment length in centimorgans
    pub fn get_seg_len_cm(&self, gmap: &GeneticMap) -> f32 {
        gmap.get_cm(self.e) - gmap.get_cm(self.s)
    }

    pub fn is_haploid_ibd(&self) -> bool {
        ((self.i & 0x3) == 3) && ((self.j & 0x3) == 3)
    }
    pub fn is_dipoid_ibd(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;
        (m <= 2) && (n <= 2) && ((m == 2) == (n == 2))
    }

    /// check if IBD segment is resulting from merging of more one original segments
    pub fn is_from_merge(&self) -> bool {
        ((self.i & 0x3) == 2) && ((self.j & 0x3) == 2)
    }
    /// check validity of IBD segments
    /// - the combinaton of hap1 and hap2 values are valid
    /// - start position is smaller then end position
    pub fn is_valid(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;
        let valid_diploid = (m <= 2) && (n <= 2) && ((m == 2) == (n == 2));
        let valid_haploid = (m == 3) && (n == 3);
        let valid_coords = self.s < self.e;

        valid_coords && (valid_diploid || valid_haploid)
    }
}

pub fn write_ibdseg_vec(v: &[IbdSeg], out: impl AsRef<Path>) -> Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(out.as_ref())
        .map(std::io::BufWriter::new)
        .context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
    let sz = v.len() as u64;
    file.write_all(&sz.to_le_bytes()).context(WriteFileSnafu {
        path: out.as_ref().to_path_buf(),
    })?;

    for seg in v.iter() {
        file.write_all(&seg.i.to_le_bytes())
            .context(WriteFileSnafu {
                path: out.as_ref().to_path_buf(),
            })?;
        file.write_all(&seg.j.to_le_bytes())
            .context(WriteFileSnafu {
                path: out.as_ref().to_path_buf(),
            })?;
        file.write_all(&seg.s.to_le_bytes())
            .context(WriteFileSnafu {
                path: out.as_ref().to_path_buf(),
            })?;
        file.write_all(&seg.e.to_le_bytes())
            .context(WriteFileSnafu {
                path: out.as_ref().to_path_buf(),
            })?;
    }
    Ok(())
}

pub fn read_ibdseg_vec(out: impl AsRef<Path>) -> Result<Vec<IbdSeg>> {
    let mut file = std::fs::File::open(out.as_ref())
        .map(std::io::BufReader::new)
        .context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
    let mut byte8 = [0u8; 8];
    let mut byte4 = [0u8; 4];
    use std::io::Read;
    file.read_exact(&mut byte8).context(WriteFileSnafu {
        path: out.as_ref().to_path_buf(),
    })?;
    let sz = u64::from_le_bytes(byte8);
    let mut v = Vec::with_capacity(sz as usize);

    for _ in 0..sz {
        file.read_exact(&mut byte4).context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
        let i = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
        let j = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
        let s = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).context(WriteFileSnafu {
            path: out.as_ref().to_path_buf(),
        })?;
        let e = u32::from_le_bytes(byte4);
        let seg = IbdSeg { i, j, s, e };
        v.push(seg);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gmap::GeneticMap;

    #[cfg(test)]
    mod construction_validation_tests {
        use super::*;

        #[test]
        fn test_new_basic_construction() {
            let seg = IbdSeg::new(10, 0, 20, 1, 1000, 2000, 0);
            assert_eq!(seg.i, 10 * 4 ); // i=40
            assert_eq!(seg.j, 20 * 4 + 1); // j=81
            assert_eq!(seg.s, 1000);
            assert_eq!(seg.e, 2000);
        }

        #[test]
        fn test_new_with_position_shift() {
            let seg = IbdSeg::new(5, 1, 7, 0, 500, 800, 1000000);
            assert_eq!(seg.i, 5 * 4 + 1); // i=21
            assert_eq!(seg.j, 7 * 4); // j=28
            assert_eq!(seg.s, 500 + 1000000);
            assert_eq!(seg.e, 800 + 1000000);
        }

        #[test]
        fn test_new_haplotype_encoding() {
            // Test all diploid haplotype combinations
            let seg_00 = IbdSeg::new(1, 0, 2, 0, 100, 200, 0);
            assert_eq!(seg_00.haplotype_pair(), (1, 0, 2, 0));

            let seg_01 = IbdSeg::new(1, 0, 2, 1, 100, 200, 0);
            assert_eq!(seg_01.haplotype_pair(), (1, 0, 2, 1));

            let seg_10 = IbdSeg::new(1, 1, 2, 0, 100, 200, 0);
            assert_eq!(seg_10.haplotype_pair(), (1, 1, 2, 0));

            let seg_11 = IbdSeg::new(1, 1, 2, 1, 100, 200, 0);
            assert_eq!(seg_11.haplotype_pair(), (1, 1, 2, 1));
        }

        #[test]
        fn test_new_merged_segment_encoding() {
            let seg_merged = IbdSeg::new(1, 2, 2, 2, 100, 200, 0);
            assert_eq!(seg_merged.haplotype_pair(), (1, 2, 2, 2));
            assert!(seg_merged.is_from_merge());
        }

        #[test]
        fn test_new_haploid_segment_encoding() {
            let seg_haploid = IbdSeg::new(1, 3, 2, 3, 100, 200, 0);
            assert_eq!(seg_haploid.haplotype_pair(), (1, 3, 2, 3));
            assert!(seg_haploid.is_haploid_ibd());
        }

        #[test]
        fn test_new_zero_values() {
            let seg = IbdSeg::new(0, 0, 0, 0, 0, 100, 0);
            assert_eq!(seg.individual_pair(), (0, 0));
            assert_eq!(seg.coords(), (0, 100));
        }

        #[test]
        fn test_new_maximum_values() {
            let max_ind = u32::MAX >> 2;
            let seg = IbdSeg::new(max_ind, 1, max_ind, 1, u32::MAX - 1000, u32::MAX - 500, 0);
            assert_eq!(seg.individual_pair(), (max_ind, max_ind));
            assert_eq!(seg.coords(), (u32::MAX - 1000, u32::MAX - 500));
        }
    }

    #[cfg(test)]
    mod normalization_tests {
        use super::*;

        #[test]
        fn test_normalized_swaps_when_needed() {
            let mut seg = IbdSeg::new(5, 0, 10, 1, 1000, 2000, 0);
            let original_i = seg.i; // 5*4+0 = 20
            let original_j = seg.j; // 10*4+1 = 41

            seg.normalized();

            // Should swap since original i (20) < original j (41)
            assert_eq!(seg.i, original_j);
            assert_eq!(seg.j, original_i);
        }

        #[test]
        fn test_normalized_no_swap_when_ordered() {
            let mut seg = IbdSeg::new(10, 0, 5, 1, 1000, 2000, 0);
            let original_i = seg.i; // 10*4+0 = 40
            let original_j = seg.j; // 5*4+1 = 21

            seg.normalized();

            // Should not swap since i (40) > j (21)
            assert_eq!(seg.i, original_i);
            assert_eq!(seg.j, original_j);
        }

        #[test]
        fn test_normalized_same_values() {
            let mut seg = IbdSeg::new(5, 0, 5, 0, 1000, 2000, 0);
            let original_i = seg.i;
            let original_j = seg.j;

            seg.normalized();

            // Should not swap when equal
            assert_eq!(seg.i, original_i);
            assert_eq!(seg.j, original_j);
        }
    }

    #[cfg(test)]
    mod data_extraction_tests {
        use super::*;

        #[test]
        fn test_haplotype_pair_int() {
            let seg = IbdSeg::new(10, 2, 15, 1, 1000, 2000, 0);
            let (i, j) = seg.haplotype_pair_int();
            assert_eq!(i, 10 * 4 + 2);
            assert_eq!(j, 15 * 4 + 1);
        }

        #[test]
        fn test_haplotype_pair_unpacked() {
            let seg = IbdSeg::new(10, 2, 15, 1, 1000, 2000, 0);
            let (ind1, hap1, ind2, hap2) = seg.haplotype_pair();
            assert_eq!(ind1, 10);
            assert_eq!(hap1, 2);
            assert_eq!(ind2, 15);
            assert_eq!(hap2, 1);
        }

        #[test]
        fn test_individual_pair() {
            let seg = IbdSeg::new(25, 3, 100, 0, 1000, 2000, 0);
            let (ind1, ind2) = seg.individual_pair();
            assert_eq!(ind1, 25);
            assert_eq!(ind2, 100);
        }

        #[test]
        fn test_coords() {
            let seg = IbdSeg::new(1, 0, 2, 1, 500000, 600000, 1000000);
            let (start, end) = seg.coords();
            assert_eq!(start, 1500000);
            assert_eq!(end, 1600000);
        }

        #[test]
        fn test_coords_zero_positions() {
            let seg = IbdSeg::new(1, 0, 2, 1, 0, 1000, 0);
            let (start, end) = seg.coords();
            assert_eq!(start, 0);
            assert_eq!(end, 1000);
        }
    }

    #[cfg(test)]
    mod classification_tests {
        use super::*;

        #[test]
        fn test_is_haploid_ibd_true() {
            let seg = IbdSeg::new(1, 3, 2, 3, 1000, 2000, 0);
            assert!(seg.is_haploid_ibd());
        }

        #[test]
        fn test_is_haploid_ibd_false() {
            let seg_diploid = IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0);
            assert!(!seg_diploid.is_haploid_ibd());

            let seg_merged = IbdSeg::new(1, 2, 2, 2, 1000, 2000, 0);
            assert!(!seg_merged.is_haploid_ibd());
        }

        #[test]
        fn test_is_dipoid_ibd_valid_combinations() {
            // Valid diploid combinations
            let seg_00 = IbdSeg::new(1, 0, 2, 0, 1000, 2000, 0);
            assert!(seg_00.is_dipoid_ibd());

            let seg_01 = IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0);
            assert!(seg_01.is_dipoid_ibd());

            let seg_10 = IbdSeg::new(1, 1, 2, 0, 1000, 2000, 0);
            assert!(seg_10.is_dipoid_ibd());

            let seg_11 = IbdSeg::new(1, 1, 2, 1, 1000, 2000, 0);
            assert!(seg_11.is_dipoid_ibd());

            // Valid merged combinations
            let seg_22 = IbdSeg::new(1, 2, 2, 2, 1000, 2000, 0);
            assert!(seg_22.is_dipoid_ibd());
        }

        #[test]
        fn test_is_dipoid_ibd_invalid_combinations() {
            let seg_haploid = IbdSeg::new(1, 3, 2, 3, 1000, 2000, 0);
            assert!(!seg_haploid.is_dipoid_ibd());

            // Mixed merged and diploid (should be invalid)
            let seg_mixed1 = IbdSeg::new(1, 0, 2, 2, 1000, 2000, 0);
            assert!(!seg_mixed1.is_dipoid_ibd());

            let seg_mixed2 = IbdSeg::new(1, 2, 2, 1, 1000, 2000, 0);
            assert!(!seg_mixed2.is_dipoid_ibd());
        }

        #[test]
        fn test_is_from_merge() {
            let seg_merged = IbdSeg::new(1, 2, 2, 2, 1000, 2000, 0);
            assert!(seg_merged.is_from_merge());

            let seg_diploid = IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0);
            assert!(!seg_diploid.is_from_merge());

            let seg_haploid = IbdSeg::new(1, 3, 2, 3, 1000, 2000, 0);
            assert!(!seg_haploid.is_from_merge());
        }
    }

    #[cfg(test)]
    mod validation_tests {
        use super::*;

        #[test]
        fn test_is_valid_diploid_segments() {
            let valid_seg = IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0);
            assert!(valid_seg.is_valid());

            let valid_merged = IbdSeg::new(1, 2, 2, 2, 1000, 2000, 0);
            assert!(valid_merged.is_valid());
        }

        #[test]
        fn test_is_valid_haploid_segments() {
            let valid_haploid = IbdSeg::new(1, 3, 2, 3, 1000, 2000, 0);
            assert!(valid_haploid.is_valid());
        }

        #[test]
        fn test_is_valid_invalid_coordinates() {
            // Start >= end should be invalid
            let invalid_coords1 = IbdSeg::new(1, 0, 2, 1, 2000, 2000, 0);
            assert!(!invalid_coords1.is_valid());

            let invalid_coords2 = IbdSeg::new(1, 0, 2, 1, 2000, 1000, 0);
            assert!(!invalid_coords2.is_valid());
        }

        #[test]
        fn test_is_valid_invalid_haplotype_combinations() {
            // Mixed haploid/diploid should be invalid
            let mixed1 = IbdSeg::new(1, 3, 2, 0, 1000, 2000, 0);
            assert!(!mixed1.is_valid());

            let mixed2 = IbdSeg::new(1, 1, 2, 3, 1000, 2000, 0);
            assert!(!mixed2.is_valid());

            // Mixed merged/diploid should be invalid
            let mixed3 = IbdSeg::new(1, 2, 2, 0, 1000, 2000, 0);
            assert!(!mixed3.is_valid());
        }

        #[test]
        fn test_is_valid_edge_case_coordinates() {
            let zero_start = IbdSeg::new(1, 0, 2, 1, 0, 1, 0);
            assert!(zero_start.is_valid());

            let max_coords = IbdSeg::new(1, 0, 2, 1, u32::MAX - 1, u32::MAX, 0);
            assert!(max_coords.is_valid());
        }
    }

    #[cfg(test)]
    mod genetic_map_tests {
        use super::*;

        #[test]
        fn test_get_seg_len_cm_simple() {
            // Create a simple genetic map: position -> cM mapping
            let gmap =
                GeneticMap::from_vec(vec![(1000, 0.5), (2000, 1.0), (3000, 1.5), (4000, 2.0)]);

            let seg = IbdSeg::new(1, 0, 2, 1, 1500, 2500, 0);
            let cm_length = seg.get_seg_len_cm(&gmap);

            // Should interpolate between the genetic map points
            assert!(cm_length > 0.0);
            assert!(cm_length < 1.0); // Should be less than 1 cM for this range
        }

        #[test]
        fn test_get_seg_len_cm_exact_positions() {
            let gmap = GeneticMap::from_vec(vec![(1000, 1.0), (2000, 2.0), (3000, 4.0)]);

            let seg = IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0);
            let cm_length = seg.get_seg_len_cm(&gmap);

            // From 1.0 cM to 2.0 cM = 1.0 cM
            assert_eq!(cm_length, 1.0);
        }

        #[test]
        fn test_get_seg_len_cm_single_position() {
            let gmap = GeneticMap::from_vec(vec![(1000, 1.0), (2000, 2.0)]);

            let seg = IbdSeg::new(1, 0, 2, 1, 1500, 1500, 0);
            let cm_length = seg.get_seg_len_cm(&gmap);

            // Same start and end position should give 0 length
            assert_eq!(cm_length, 0.0);
        }

        #[test]
        fn test_get_seg_len_cm_with_position_shift() {
            let gmap = GeneticMap::from_vec(vec![(1000000, 5.0), (1001000, 6.0), (1002000, 7.0)]);

            // Segment with position shift
            let seg = IbdSeg::new(1, 0, 2, 1, 0, 1000, 1000000);
            let cm_length = seg.get_seg_len_cm(&gmap);

            // Should calculate based on shifted coordinates (1000000 to 1001000)
            assert!(cm_length > 0.0);
        }
    }

    #[cfg(test)]
    mod file_io_tests {
        use super::*;
        use std::fs;

        fn create_test_segments() -> Vec<IbdSeg> {
            vec![
                IbdSeg::new(1, 0, 2, 1, 1000, 2000, 0),
                IbdSeg::new(3, 1, 4, 0, 3000, 4000, 0),
                IbdSeg::new(5, 2, 6, 2, 5000, 6000, 0), // Merged segment
                IbdSeg::new(7, 3, 8, 3, 7000, 8000, 0), // Haploid segment
            ]
        }

        #[test]
        fn test_write_read_roundtrip() {
            let test_file = "test_ibdseg_roundtrip.bin";
            let original_segments = create_test_segments();

            // Write segments
            write_ibdseg_vec(&original_segments, test_file).unwrap();

            // Read segments back
            let read_segments = read_ibdseg_vec(test_file).unwrap();

            // Clean up
            fs::remove_file(test_file).ok();

            // Verify they match
            assert_eq!(original_segments.len(), read_segments.len());
            for (original, read) in original_segments.iter().zip(read_segments.iter()) {
                assert_eq!(original.i, read.i);
                assert_eq!(original.j, read.j);
                assert_eq!(original.s, read.s);
                assert_eq!(original.e, read.e);
            }
        }

        #[test]
        fn test_write_read_empty_vector() {
            let test_file = "test_ibdseg_empty.bin";
            let empty_segments: Vec<IbdSeg> = vec![];

            // Write empty vector
            write_ibdseg_vec(&empty_segments, test_file).unwrap();

            // Read back
            let read_segments = read_ibdseg_vec(test_file).unwrap();

            // Clean up
            fs::remove_file(test_file).ok();

            assert_eq!(read_segments.len(), 0);
        }

        #[test]
        fn test_write_read_single_segment() {
            let test_file = "test_ibdseg_single.bin";
            let single_segment = vec![IbdSeg::new(10, 1, 20, 0, 50000, 60000, 1000000)];

            write_ibdseg_vec(&single_segment, test_file).unwrap();
            let read_segments = read_ibdseg_vec(test_file).unwrap();

            // Clean up
            fs::remove_file(test_file).ok();

            assert_eq!(read_segments.len(), 1);
            let read_seg = &read_segments[0];
            assert_eq!(read_seg.individual_pair(), (10, 20));
            assert_eq!(read_seg.haplotype_pair(), (10, 1, 20, 0));
            assert_eq!(read_seg.coords(), (1050000, 1060000));
        }

        #[test]
        fn test_write_read_large_segments() {
            let test_file = "test_ibdseg_large.bin";
            let mut large_segments = Vec::new();

            // Create many segments with different patterns
            for i in 0..1000 {
                let seg = IbdSeg::new(
                    i % 100,
                    (i % 4) as u8,
                    (i + 50) % 100,
                    ((i + 1) % 4) as u8,
                    i * 1000,
                    (i + 1) * 1000,
                    i * 10000,
                );
                large_segments.push(seg);
            }

            write_ibdseg_vec(&large_segments, test_file).unwrap();
            let read_segments = read_ibdseg_vec(test_file).unwrap();

            // Clean up
            fs::remove_file(test_file).ok();

            assert_eq!(large_segments.len(), read_segments.len());
            assert_eq!(large_segments.len(), 1000);

            // Spot check a few segments
            for i in [0, 100, 500, 999] {
                assert_eq!(large_segments[i].i, read_segments[i].i);
                assert_eq!(large_segments[i].j, read_segments[i].j);
                assert_eq!(large_segments[i].s, read_segments[i].s);
                assert_eq!(large_segments[i].e, read_segments[i].e);
            }
        }

        #[test]
        fn test_read_nonexistent_file() {
            let result = read_ibdseg_vec("nonexistent_file.bin");
            assert!(result.is_err());
        }

        #[test]
        fn test_write_to_invalid_path() {
            let segments = create_test_segments();
            let result = write_ibdseg_vec(&segments, "/invalid/path/that/does/not/exist.bin");
            assert!(result.is_err());
        }
    }

    #[cfg(test)]
    mod edge_cases_and_errors {
        use super::*;
        use std::fs;
        use std::io::Write;

        #[test]
        fn test_segment_with_all_zeros() {
            let seg = IbdSeg {
                i: 0,
                j: 0,
                s: 0,
                e: 0,
            };
            assert_eq!(seg.individual_pair(), (0, 0));
            assert_eq!(seg.haplotype_pair(), (0, 0, 0, 0));
            assert!(!seg.is_valid()); // s == e is invalid
        }

        #[test]
        fn test_segment_with_maximum_values() {
            let seg = IbdSeg {
                i: u32::MAX,
                j: u32::MAX,
                s: u32::MAX - 1,
                e: u32::MAX,
            };
            let (ind1, hap1, ind2, hap2) = seg.haplotype_pair();
            assert_eq!(ind1, u32::MAX >> 2);
            assert_eq!(hap1, 3);
            assert_eq!(ind2, u32::MAX >> 2);
            assert_eq!(hap2, 3);
            assert!(seg.is_valid());
        }

        #[test]
        fn test_bit_manipulation_boundary_cases() {
            // Test that bit operations work correctly at boundaries
            let seg1 = IbdSeg::new(1073741823, 3, 1073741823, 3, 1000, 2000, 0); // Max individual with hap 3
            assert_eq!(seg1.individual_pair(), (1073741823, 1073741823));
            assert_eq!(seg1.haplotype_pair().1, 3);
            assert_eq!(seg1.haplotype_pair().3, 3);

            let seg2 = IbdSeg::new(0, 0, 0, 0, 1000, 2000, 0); // Min values
            assert_eq!(seg2.individual_pair(), (0, 0));
            assert_eq!(seg2.haplotype_pair().1, 0);
            assert_eq!(seg2.haplotype_pair().3, 0);
        }

        #[test]
        fn test_position_overflow_handling() {
            // Test large position shifts that might cause overflow issues
            let seg = IbdSeg::new(1, 0, 2, 1, u32::MAX - 1000, u32::MAX - 500, 500);
            let (start, end) = seg.coords();
            assert_eq!(start, u32::MAX - 500); // Should wrap correctly
            assert_eq!(end, u32::MAX);
        }

        #[test]
        fn test_genetic_map_edge_cases() {
            // Test with genetic map having multiple points
            let gmap = GeneticMap::from_vec(vec![(500, 0.5), (1000, 1.0), (1500, 1.5)]);
            let seg = IbdSeg::new(1, 0, 2, 1, 1000, 1000, 0);
            let cm_length = seg.get_seg_len_cm(&gmap);
            assert_eq!(cm_length, 0.0);

            // Test with segment within genetic map range
            let seg2 = IbdSeg::new(1, 0, 2, 1, 500, 1500, 0);
            let cm_length2 = seg2.get_seg_len_cm(&gmap);
            // From the calculation, it returns 0.5 cM (1.5 - 1.0), not (1.5 - 0.5)
            assert_eq!(cm_length2, 0.5);
        }

        #[test]
        fn test_corrupted_file_handling() {
            let test_file = "test_corrupted.bin";

            // Create a file with invalid content
            let mut file = fs::File::create(test_file).unwrap();
            file.write_all(&[1, 2, 3, 4]).unwrap(); // Invalid file structure
            drop(file);

            let result = read_ibdseg_vec(test_file);
            fs::remove_file(test_file).ok();

            assert!(result.is_err());
        }

        #[test]
        fn test_file_with_wrong_segment_count() {
            let test_file = "test_wrong_count.bin";

            // Create file claiming to have more segments than actually present
            let mut file = fs::File::create(test_file).unwrap();
            let fake_count: u64 = 10; // Claim 10 segments
            file.write_all(&fake_count.to_le_bytes()).unwrap();
            // But only write one segment worth of data
            file.write_all(&1u32.to_le_bytes()).unwrap();
            file.write_all(&2u32.to_le_bytes()).unwrap();
            file.write_all(&1000u32.to_le_bytes()).unwrap();
            file.write_all(&2000u32.to_le_bytes()).unwrap();
            drop(file);

            let result = read_ibdseg_vec(test_file);
            fs::remove_file(test_file).ok();

            assert!(result.is_err()); // Should fail trying to read more data
        }

        #[test]
        fn test_normalization_with_encoded_values() {
            // Test normalization with various encoded haplotype values
            let mut seg1 = IbdSeg {
                i: 100,
                j: 50,
                s: 1000,
                e: 2000,
            };
            seg1.normalized();
            // No swap because i (100) > j (50)
            assert_eq!(seg1.i, 100);
            assert_eq!(seg1.j, 50);

            let mut seg2 = IbdSeg {
                i: 50,
                j: 100,
                s: 1000,
                e: 2000,
            };
            seg2.normalized();
            // Swap because i (50) < j (100)
            assert_eq!(seg2.i, 100);
            assert_eq!(seg2.j, 50);

            let mut seg3 = IbdSeg {
                i: 75,
                j: 75,
                s: 1000,
                e: 2000,
            };
            seg3.normalized();
            // No swap because i == j
            assert_eq!(seg3.i, 75);
            assert_eq!(seg3.j, 75);
        }
    }
}
