use arrow_array::builder::BooleanBufferBuilder;
use arrow_array::ArrayRef;
use arrow_array::BooleanArray;
use arrow_array::RecordBatch;
use arrow_buffer::BooleanBuffer;
use arrow_schema::ArrowError;

use snafu::prelude::*;

use crate::site::Sites;
use bitvec::prelude::*;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};
use std::backtrace::Backtrace;
use std::fs::File;
use std::num::ParseIntError;
use std::path::Path;
use std::sync::Arc;

// Matrix of 0 and 1.
///
/// Multiallelic sites must be presented as multiple lines of biallelic genotypes
#[derive(Debug, Serialize, Deserialize)]
pub struct GenotypeMatrix {
    bv: BitVec<u64, Lsb0>,
    ncols: usize,
}

type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Snafu)]
pub enum Error {
    Arrow {
        #[snafu(source(from(ArrowError, Box::new)))]
        source: Box<ArrowError>,
        backtrace: Box<Option<Backtrace>>,
    },
    Parquet {
        #[snafu(source(from(parquet::errors::ParquetError, Box::new)))]
        source: Box<parquet::errors::ParquetError>,
        backtrace: Box<Option<Backtrace>>,
    },
    Io {
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    ParseInt {
        source: ParseIntError,
        backtrace: Box<Option<Backtrace>>,
    },
    DowncastToBoolArray {
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Row order length {row_order_len} does not match matrix rows {matrix_rows}"))]
    RowOrderLengthMismatch {
        row_order_len: usize,
        matrix_rows: usize,
        backtrace: Box<Option<Backtrace>>,
    },
    MissingGenotype {
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("BitVec length mismatch after append: expected {expected}, got {actual}"))]
    BitVecLengthMismatch {
        expected: usize,
        actual: usize,
        backtrace: Box<Option<Backtrace>>,
    },
}

impl GenotypeMatrix {
    pub fn new(ncol: usize) -> Self {
        Self {
            bv: BitVec::<u64, Lsb0>::new(),
            ncols: ncol,
        }
    }

    pub fn nrows(&self) -> usize {
        if self.ncols == 0 {
            0
        } else {
            let nrow = self.bv.len() / self.ncols;
            assert_eq!(self.ncols * nrow, self.bv.len());
            nrow
        }
    }

    pub fn ncols(&self) -> usize {
        self.ncols
    }

    pub fn extend_gt_calls(&mut self, gt_bool_it: impl Iterator<Item = bool>) {
        self.bv.extend(gt_bool_it);
    }

    pub fn get_at(&self, pos_idx: usize, genome_id: usize) -> bool {
        let index = pos_idx * self.ncols + genome_id;
        self.bv[index]
    }

    pub fn get_row(&self, pos_idx: usize) -> &BitSlice<u64> {
        let s = pos_idx * self.ncols;
        let e = s + self.ncols;
        &self.bv[s..e]
    }

    pub fn get_slice(&self, row: usize, col1: usize, col2: usize) -> &BitSlice<u64, Lsb0> {
        &self.bv[(row * self.ncols + col1)..(row * self.ncols + col2)]
    }

    pub fn has_too_many_discod_sites(
        &self,
        ind_pair: (u32, u32),
        col1: usize,
        col2: usize,
        max_ndiscord: u32,
    ) -> bool {
        let s1 = self.get_slice(ind_pair.0 as usize * 2, col1, col2);
        let s2 = self.get_slice(ind_pair.0 as usize * 2 + 1, col1, col2);
        let s3 = self.get_slice(ind_pair.1 as usize * 2, col1, col2);
        let s4 = self.get_slice(ind_pair.1 as usize * 2 + 1, col1, col2);

        let mut ndicosrd = 0u32;
        for (a, b, c, d) in itertools::multizip((s1, s2, s3, s4)) {
            let a = *a;
            let b = *b;
            let c = *c;
            let d = *d;
            // if (a != b) && (a != c) && (a != d) && (b != c) && (b != d) && (c != d) {
            if (a == b) && (c == d) && (a != c) {
                ndicosrd += 1;
                if ndicosrd > max_ndiscord {
                    return true;
                }
            }
        }

        false
    }

    pub fn transpose(&self) -> Self {
        let mut bv2 = BitVec::with_capacity(self.bv.len());
        for col in 0..self.ncols() {
            for row in 0..self.nrows() {
                bv2.push(self.get_at(row, col))
            }
        }
        let ncols2 = self.nrows();
        Self {
            bv: bv2,
            ncols: ncols2,
        }
    }

    pub fn append_row(&mut self, other: &BitSlice<u64>) -> Result<()> {
        let len1 = self.bv.len();
        self.bv.extend_from_bitslice(other);
        let len2 = self.bv.len();
        ensure!(
            len1 + self.ncols == len2,
            BitVecLengthMismatchSnafu {
                expected: len1 + self.ncols,
                actual: len2
            }
        );
        Ok(())
    }

    pub fn reorder_rows(self, row_orders: &[u32]) -> Result<Self> {
        ensure!(
            row_orders.len() == self.nrows(),
            RowOrderLengthMismatchSnafu {
                row_order_len: row_orders.len(),
                matrix_rows: self.nrows()
            }
        );
        let mut gt2 = GenotypeMatrix::new(self.ncols);
        gt2.bv.reserve(self.bv.len());

        for idx in row_orders {
            let row_slice = self.get_row(*idx as usize);
            gt2.append_row(row_slice)?;
        }
        Ok(gt2)
    }

    pub fn merge(&mut self, other: Self) {
        assert_eq!(self.ncols, other.ncols);
        let slice = other.bv.as_bitslice();
        self.bv.extend_from_bitslice(slice);
    }

    /// Checks if a genotype contains a specific allele at a given genomic position.
    ///
    /// This method determines whether an individual (genotype ID) carries a particular
    /// allele at a specified genomic position. It handles both reference alleles
    /// (encoded as "REF") and alternative alleles.
    ///
    /// # Arguments
    ///
    /// * `pos` - Genomic position (in genome-wide coordinates)
    /// * `gid` - Genotype ID (individual identifier)
    /// * `allele` - Allele sequence as bytes (use b"REF" for reference allele)
    /// * `sites` - Sites container with variant information and allele mappings
    ///
    /// # Returns
    ///
    /// `true` if the individual carries the specified allele at the position,
    /// `false` otherwise
    ///
    /// # Behavior
    ///
    /// - For reference alleles (b"REF"): Returns `true` if ALL sites in the position
    ///   range are homozygous reference (no alternative alleles present)
    /// - For alternative alleles: Returns `true` if ANY site in the position range
    ///   matches the allele sequence AND the individual carries that allele
    /// - Returns `false` if the position is not found in the sites container
    pub fn has_allele(&self, pos: u32, gid: u32, allele: &[u8], sites: &Sites) -> bool {
        let (s, e) = sites.get_idx_by_position(pos);
        if s == e {
            false
        } else {
            match allele {
                b"REF" => (s..e).all(|pos_idx| !self.get_at(pos_idx, gid as usize)),
                _ => (s..e).any(|pos_idx| {
                    let stored_alleles = sites.get_alleles_by_idx(pos_idx);
                    // Check if the target allele exists in the space-separated alleles
                    let has_target_allele = stored_alleles
                        .split(|&x| x == b' ')
                        .any(|stored_allele| stored_allele == allele);
                    has_target_allele && self.get_at(pos_idx, gid as usize)
                }),
            }
        }
    }
    pub fn into_parquet_file(self, p: impl AsRef<Path>) -> Result<()> {
        // build array from genotype matrix
        let mut builder = BooleanBufferBuilder::new(self.bv.len());
        for i in 0..self.bv.len() {
            builder.append(self.bv[i]);
        }

        let buffer: BooleanBuffer = builder.finish();
        let arr = BooleanArray::new(buffer, None);

        // trick: use field name to save the ncol member variable
        let fieldname = format!("{}", self.ncols);

        let batch = RecordBatch::try_from_iter(vec![(fieldname, Arc::new(arr) as ArrayRef)])
            .context(ArrowSnafu {})?;
        // writer
        let file = File::create(p.as_ref()).context(IoSnafu {})?;
        // -- default writer properties
        let props = WriterProperties::builder().build();
        let mut writer =
            ArrowWriter::try_new(file, batch.schema(), Some(props)).context(ParquetSnafu {})?;
        // write batch
        writer.write(&batch).context(ParquetSnafu {})?;
        // writer must be closed to write footer
        writer.close().context(ParquetSnafu {})?;
        Ok(())
    }

    pub fn from_parquet_file(p: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(p).context(IoSnafu {})?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).context(ParquetSnafu {})?;
        let ncol = builder.schema().fields[0]
            .name()
            .parse::<usize>()
            .context(ParseIntSnafu {})?;
        let mut reader = builder.build().context(ParquetSnafu {})?;
        let mut gm = GenotypeMatrix::new(ncol);
        for record_batch in &mut reader {
            let record_batch = record_batch.context(ArrowSnafu {})?;
            record_batch
                .column(0)
                .as_any()
                .downcast_ref::<BooleanArray>()
                .context(DowncastToBoolArraySnafu {})?
                .into_iter()
                .try_for_each(|x| match x {
                    Some(x) => {
                        gm.bv.push(x);
                        Ok(())
                    }
                    _ => MissingGenotypeSnafu {}.fail(),
                })?;
        }

        Ok(gm)
    }

    /// calculate allele frequency assume each row of the matrix is a biallelic site, and each column is a haplotype
    pub fn get_afreq(&self) -> Vec<f64> {
        let mut afreq = Vec::<f64>::with_capacity(self.nrows());
        for r in 0..self.nrows() {
            let ac = self.get_row(r).iter().map(|x| *x as usize).sum::<usize>() as f64;
            let n = self.ncols() as f64;
            afreq.push(ac / n);
        }
        afreq
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_gt_matrix() {
        // rows are genomes and cols are sites
        let mut m = GenotypeMatrix::new(5);
        let gt = [
            [1, 0, 1, 0, 1],
            [1, 0, 1, 0, 1],
            [1, 0, 0, 1, 1],
            [1, 0, 0, 1, 1],
        ];
        let it = gt.into_iter().flat_map(|s| s.map(|x| x == 1).into_iter());
        m.extend_gt_calls(it);

        // two many
        assert!(m.has_too_many_discod_sites((0, 1), 0, 5, 1));
        // not too many
        assert!(!m.has_too_many_discod_sites((0, 1), 0, 5, 2));
        assert!(!m.has_too_many_discod_sites((0, 1), 0, 5, 3));
    }

    mod construction {
        use super::*;

        #[test]
        fn test_new_empty_matrix() {
            let m = GenotypeMatrix::new(10);
            assert_eq!(m.ncols(), 10);
            assert_eq!(m.nrows(), 0);
            assert_eq!(m.bv.len(), 0);
        }

        #[test]
        fn test_extend_gt_calls_single_row() {
            let mut m = GenotypeMatrix::new(3);
            m.extend_gt_calls([true, false, true].into_iter());
            assert_eq!(m.nrows(), 1);
            assert_eq!(m.ncols(), 3);
            assert!(m.get_at(0, 0));
            assert!(!m.get_at(0, 1));
            assert!(m.get_at(0, 2));
        }

        #[test]
        fn test_extend_gt_calls_multiple_rows() {
            let mut m = GenotypeMatrix::new(2);
            m.extend_gt_calls([true, false, false, true].into_iter());
            assert_eq!(m.nrows(), 2);
            assert_eq!(m.ncols(), 2);
            assert!(m.get_at(0, 0));
            assert!(!m.get_at(0, 1));
            assert!(!m.get_at(1, 0));
            assert!(m.get_at(1, 1));
        }
    }

    mod access_methods {
        use super::*;

        fn create_test_matrix() -> GenotypeMatrix {
            let mut m = GenotypeMatrix::new(3);
            m.extend_gt_calls(
                [
                    true, false, true, // Row 0
                    false, true, false, // Row 1
                    true, true, true, // Row 2
                ]
                .into_iter(),
            );
            m
        }

        #[test]
        fn test_get_at() {
            let m = create_test_matrix();
            assert!(m.get_at(0, 0));
            assert!(!m.get_at(0, 1));
            assert!(m.get_at(0, 2));
            assert!(!m.get_at(1, 0));
            assert!(m.get_at(1, 1));
            assert!(!m.get_at(1, 2));
            assert!(m.get_at(2, 0));
            assert!(m.get_at(2, 1));
            assert!(m.get_at(2, 2));
        }

        #[test]
        fn test_get_row() {
            let m = create_test_matrix();
            let row0 = m.get_row(0);
            assert_eq!(row0.len(), 3);
            assert!(row0[0]);
            assert!(!row0[1]);
            assert!(row0[2]);

            let row1 = m.get_row(1);
            assert!(!row1[0]);
            assert!(row1[1]);
            assert!(!row1[2]);
        }

        #[test]
        fn test_get_slice() {
            let m = create_test_matrix();
            let slice = m.get_slice(0, 1, 3);
            assert_eq!(slice.len(), 2);
            assert!(!slice[0]); // m.get_at(0, 1)
            assert!(slice[1]); // m.get_at(0, 2)
        }

        #[test]
        #[should_panic]
        fn test_get_at_out_of_bounds() {
            let m = create_test_matrix();
            m.get_at(3, 0); // Row 3 doesn't exist
        }

        #[test]
        fn test_get_at_col_bounds() {
            let m = create_test_matrix(); // This has ncols=3, so valid indices are 0,1,2
                                          // Valid accesses should work
            assert!(m.get_at(0, 0));
            assert!(!m.get_at(0, 1));
            assert!(m.get_at(0, 2));

            // Out of bounds access should panic - test with clearly out of bounds index
            let result = std::panic::catch_unwind(|| {
                m.get_at(0, 10); // This should definitely panic
            });
            assert!(result.is_err());
        }
    }

    mod matrix_operations {
        use super::*;

        #[test]
        fn test_transpose() {
            let mut m = GenotypeMatrix::new(3);
            m.extend_gt_calls(
                [
                    true, false, true, // Row 0
                    false, true, false, // Row 1
                ]
                .into_iter(),
            );

            let t = m.transpose();
            assert_eq!(t.nrows(), 3);
            assert_eq!(t.ncols(), 2);

            // Original (0,0) -> Transposed (0,0)
            assert!(t.get_at(0, 0));
            // Original (0,1) -> Transposed (1,0)
            assert!(!t.get_at(1, 0));
            // Original (0,2) -> Transposed (2,0)
            assert!(t.get_at(2, 0));
            // Original (1,0) -> Transposed (0,1)
            assert!(!t.get_at(0, 1));
            // Original (1,1) -> Transposed (1,1)
            assert!(t.get_at(1, 1));
            // Original (1,2) -> Transposed (2,1)
            assert!(!t.get_at(2, 1));
        }

        #[test]
        fn test_transpose_empty_matrix() {
            let m = GenotypeMatrix::new(5);
            let t = m.transpose();
            // When original matrix has 0 rows and 5 cols, transpose has 5 rows and 0 cols
            assert_eq!(t.ncols(), 0); // Original nrows becomes new ncols
            assert_eq!(t.bv.len(), 0); // Should be empty
        }

        #[test]
        fn test_append_row() {
            let mut m = GenotypeMatrix::new(3);
            m.extend_gt_calls([true, false, true].into_iter());

            let mut new_row_bits = BitVec::<u64, bitvec::order::Lsb0>::new();
            new_row_bits.push(false);
            new_row_bits.push(true);
            new_row_bits.push(false);
            m.append_row(new_row_bits.as_bitslice()).unwrap();

            assert_eq!(m.nrows(), 2);
            assert!(!m.get_at(1, 0));
            assert!(m.get_at(1, 1));
            assert!(!m.get_at(1, 2));
        }

        #[test]
        fn test_append_row_wrong_size() {
            let mut m = GenotypeMatrix::new(3);
            let mut wrong_size_row = BitVec::<u64, bitvec::order::Lsb0>::new();
            wrong_size_row.push(true);
            wrong_size_row.push(false); // Only 2 elements
            let result = m.append_row(wrong_size_row.as_bitslice());
            assert!(result.is_err());
        }

        #[test]
        fn test_reorder_rows() {
            let mut m = GenotypeMatrix::new(2);
            m.extend_gt_calls(
                [
                    true, false, // Row 0
                    false, true, // Row 1
                    true, true, // Row 2
                ]
                .into_iter(),
            );

            let reordered = m.reorder_rows(&[2, 0, 1]).unwrap();
            assert_eq!(reordered.nrows(), 3);
            assert_eq!(reordered.ncols(), 2);

            // Row 0 should now be original row 2
            assert!(reordered.get_at(0, 0));
            assert!(reordered.get_at(0, 1));

            // Row 1 should now be original row 0
            assert!(reordered.get_at(1, 0));
            assert!(!reordered.get_at(1, 1));

            // Row 2 should now be original row 1
            assert!(!reordered.get_at(2, 0));
            assert!(reordered.get_at(2, 1));
        }

        #[test]
        fn test_reorder_rows_wrong_length() {
            let mut m = GenotypeMatrix::new(2);
            m.extend_gt_calls([true, false].into_iter());
            let result = m.reorder_rows(&[0, 1]); // Matrix only has 1 row, not 2
            assert!(result.is_err());
        }

        #[test]
        fn test_merge_matrices() {
            let mut m1 = GenotypeMatrix::new(2);
            m1.extend_gt_calls([true, false, false, true].into_iter());

            let mut m2 = GenotypeMatrix::new(2);
            m2.extend_gt_calls([true, true].into_iter());

            m1.merge(m2);
            assert_eq!(m1.nrows(), 3);
            assert_eq!(m1.ncols(), 2);

            // First two rows unchanged
            assert!(m1.get_at(0, 0));
            assert!(!m1.get_at(0, 1));
            assert!(!m1.get_at(1, 0));
            assert!(m1.get_at(1, 1));

            // Third row from merged matrix
            assert!(m1.get_at(2, 0));
            assert!(m1.get_at(2, 1));
        }

        #[test]
        #[should_panic]
        fn test_merge_different_ncols() {
            let mut m1 = GenotypeMatrix::new(2);
            let m2 = GenotypeMatrix::new(3);
            m1.merge(m2);
        }
    }

    mod discord_analysis {
        use super::*;

        #[test]
        fn test_has_too_many_discod_sites_no_discord() {
            let mut m = GenotypeMatrix::new(4);
            // Create 4 haplotypes representing 2 diploid individuals
            // Rows 0,1 = individual 0's haplotypes; rows 2,3 = individual 1's haplotypes
            // Create concordant sites where both individuals have same homozygous genotype
            m.extend_gt_calls(
                [
                    true, true, true, true, // Haplotype 0 (ind 0, copy 1) - all 1s
                    true, true, true, true, // Haplotype 1 (ind 0, copy 2) - all 1s
                    true, true, true, true, // Haplotype 2 (ind 1, copy 1) - all 1s
                    true, true, true, true, // Haplotype 3 (ind 1, copy 2) - all 1s
                ]
                .into_iter(),
            );

            // Both individuals are homozygous 1,1 at all sites - no discord
            assert!(!m.has_too_many_discod_sites((0, 1), 0, 4, 0));
            assert!(!m.has_too_many_discod_sites((0, 1), 0, 4, 1));
        }

        #[test]
        fn test_has_too_many_discod_sites_with_discord() {
            let mut m = GenotypeMatrix::new(4);
            // Create 4 haplotypes representing 2 diploid individuals
            // Create discordant sites: ind0 homozygous 1,1 vs ind1 homozygous 0,0
            m.extend_gt_calls(
                [
                    true, true, true, true, // Haplotype 0 (ind 0, copy 1) - all 1s
                    true, true, true, true, // Haplotype 1 (ind 0, copy 2) - all 1s
                    false, false, false, false, // Haplotype 2 (ind 1, copy 1) - all 0s
                    false, false, false, false, // Haplotype 3 (ind 1, copy 2) - all 0s
                ]
                .into_iter(),
            );

            // All 4 sites are discordant: ind0=(1,1) vs ind1=(0,0)
            assert!(m.has_too_many_discod_sites((0, 1), 0, 4, 3)); // More than 3 discord
            assert!(!m.has_too_many_discod_sites((0, 1), 0, 4, 4)); // Allow up to 4 discord
        }

        #[test]
        fn test_has_too_many_discod_sites_partial_range() {
            let mut m = GenotypeMatrix::new(4);
            m.extend_gt_calls(
                [
                    true, true, false, false, // Haplotype 0 (ind 0, copy 1)
                    true, true, false, false, // Haplotype 1 (ind 0, copy 2)
                    false, false, true, true, // Haplotype 2 (ind 1, copy 1)
                    false, false, true, true, // Haplotype 3 (ind 1, copy 2)
                ]
                .into_iter(),
            );

            // Check only first 2 columns (should have discords at positions 0,1)
            assert!(m.has_too_many_discod_sites((0, 1), 0, 2, 1));
            // Check columns 2-4 (should also have discords)
            assert!(m.has_too_many_discod_sites((0, 1), 2, 4, 1));
        }
    }

    mod allele_frequency {
        use super::*;

        #[test]
        fn test_get_afreq_basic() {
            let mut m = GenotypeMatrix::new(4);
            m.extend_gt_calls(
                [
                    true, false, true, false, // 2/4 = 0.5
                    false, false, false, false, // 0/4 = 0.0
                    true, true, true, true, // 4/4 = 1.0
                ]
                .into_iter(),
            );

            let afreq = m.get_afreq();
            assert_eq!(afreq.len(), 3);
            assert!((afreq[0] - 0.5).abs() < 1e-10);
            assert!((afreq[1] - 0.0).abs() < 1e-10);
            assert!((afreq[2] - 1.0).abs() < 1e-10);
        }

        #[test]
        fn test_get_afreq_empty() {
            let m = GenotypeMatrix::new(5);
            let afreq = m.get_afreq();
            assert_eq!(afreq.len(), 0);
        }

        #[test]
        fn test_get_afreq_single_column() {
            let mut m = GenotypeMatrix::new(1);
            m.extend_gt_calls([true, false, true].into_iter());

            let afreq = m.get_afreq();
            assert_eq!(afreq.len(), 3);
            assert!((afreq[0] - 1.0).abs() < 1e-10);
            assert!((afreq[1] - 0.0).abs() < 1e-10);
            assert!((afreq[2] - 1.0).abs() < 1e-10);
        }
    }

    mod parquet_serialization {
        use super::*;

        #[test]
        fn test_parquet_round_trip() -> Result<()> {
            let mut m = GenotypeMatrix::new(3);
            m.extend_gt_calls([true, false, true, false, true, false].into_iter());

            let temp_file = NamedTempFile::new().expect("Failed to create temp file");
            let temp_path = temp_file.path();

            m.into_parquet_file(temp_path)?;
            let loaded_m = GenotypeMatrix::from_parquet_file(temp_path)?;

            assert_eq!(loaded_m.nrows(), 2);
            assert_eq!(loaded_m.ncols(), 3);
            assert!(loaded_m.get_at(0, 0));
            assert!(!loaded_m.get_at(0, 1));
            assert!(loaded_m.get_at(0, 2));
            assert!(!loaded_m.get_at(1, 0));
            assert!(loaded_m.get_at(1, 1));
            assert!(!loaded_m.get_at(1, 2));

            Ok(())
        }

        #[test]
        fn test_parquet_empty_matrix() -> Result<()> {
            let m = GenotypeMatrix::new(5);

            let temp_file = NamedTempFile::new().expect("Failed to create temp file");
            let temp_path = temp_file.path();

            m.into_parquet_file(temp_path)?;
            let loaded_m = GenotypeMatrix::from_parquet_file(temp_path)?;

            assert_eq!(loaded_m.nrows(), 0);
            assert_eq!(loaded_m.ncols(), 5);

            Ok(())
        }

        #[test]
        fn test_parquet_large_matrix() -> Result<()> {
            let mut m = GenotypeMatrix::new(100);
            // Create a pattern: alternating true/false
            let data: Vec<bool> = (0..1000).map(|i| i % 2 == 0).collect();
            m.extend_gt_calls(data.into_iter());

            let temp_file = NamedTempFile::new().expect("Failed to create temp file");
            let temp_path = temp_file.path();

            m.into_parquet_file(temp_path)?;
            let loaded_m = GenotypeMatrix::from_parquet_file(temp_path)?;

            assert_eq!(loaded_m.nrows(), 10);
            assert_eq!(loaded_m.ncols(), 100);

            // Verify pattern
            for row in 0..10 {
                for col in 0..100 {
                    let expected = (row * 100 + col) % 2 == 0;
                    assert_eq!(loaded_m.get_at(row, col), expected);
                }
            }

            Ok(())
        }
    }

    mod has_allele {
        use super::*;
        use crate::site::{AlleleBuffer, Sites};

        fn create_test_sites() -> Sites {
            let mut sites = Sites::default();
            let mut ab = AlleleBuffer::new();
            ab.push(b"REF");
            ab.push(b"T");
            sites.add_site(100, &ab).unwrap();

            sites
        }

        #[test]
        fn test_has_allele_basic() {
            let mut m = GenotypeMatrix::new(2);
            m.extend_gt_calls([true, false, false, true].into_iter());

            let sites = create_test_sites();

            assert!(m.has_allele(100, 0, b"T", &sites));
        }
    }

    mod edge_cases {
        use super::*;

        #[test]
        fn test_zero_columns() {
            let m = GenotypeMatrix::new(0);
            assert_eq!(m.ncols(), 0);
            // nrows() will be 0 when bv.len() is 0, regardless of ncols
            // The division by zero issue occurs when there are bits but ncols=0
            assert_eq!(m.bv.len(), 0);
        }

        #[test]
        fn test_single_column_single_row() {
            let mut m = GenotypeMatrix::new(1);
            m.extend_gt_calls([true].into_iter());

            assert_eq!(m.nrows(), 1);
            assert_eq!(m.ncols(), 1);
            assert!(m.get_at(0, 0));

            let row = m.get_row(0);
            assert_eq!(row.len(), 1);
            assert!(row[0]);
        }

        #[test]
        fn test_large_dimensions() {
            let ncols = 1000;
            let nrows = 500;
            let mut m = GenotypeMatrix::new(ncols);

            // Fill with alternating pattern
            let data: Vec<bool> = (0..(ncols * nrows)).map(|i| i % 2 == 0).collect();
            m.extend_gt_calls(data.into_iter());

            assert_eq!(m.nrows(), nrows);
            assert_eq!(m.ncols(), ncols);

            // Test some random positions
            assert!(m.get_at(0, 0)); // Even position
            assert!(!m.get_at(0, 1)); // Odd position
            assert!(m.get_at(1, 0)); // Even position (1000)
            assert!(!m.get_at(1, 1)); // Odd position (1001)
        }

        #[test]
        fn test_transpose_single_row() {
            let mut m = GenotypeMatrix::new(5);
            m.extend_gt_calls([true, false, true, false, true].into_iter());

            let t = m.transpose();
            assert_eq!(t.ncols(), 1); // Original nrows (1) becomes new ncols
            assert_eq!(t.nrows(), 5); // Original ncols (5) becomes new nrows

            assert!(t.get_at(0, 0));
            assert!(!t.get_at(1, 0));
            assert!(t.get_at(2, 0));
            assert!(!t.get_at(3, 0));
            assert!(t.get_at(4, 0));
        }

        #[test]
        fn test_transpose_single_column() {
            let mut m = GenotypeMatrix::new(1);
            m.extend_gt_calls([true, false, true].into_iter());

            let t = m.transpose();
            assert_eq!(t.nrows(), 1);
            assert_eq!(t.ncols(), 3);

            assert!(t.get_at(0, 0));
            assert!(!t.get_at(0, 1));
            assert!(t.get_at(0, 2));
        }
    }

    mod error_handling {
        use super::*;
        use std::fs;

        #[test]
        fn test_parquet_invalid_file() {
            let result = GenotypeMatrix::from_parquet_file("/nonexistent/path/file.parquet");
            assert!(result.is_err());
            assert!(matches!(result.unwrap_err(), Error::Io { .. }));
        }

        #[test]
        fn test_parquet_corrupted_file() {
            let temp_file = NamedTempFile::new().expect("Failed to create temp file");
            let temp_path = temp_file.path();

            // Write invalid content
            fs::write(temp_path, b"invalid parquet data").expect("Failed to write file");

            let result = GenotypeMatrix::from_parquet_file(temp_path);
            assert!(result.is_err());
            assert!(matches!(result.unwrap_err(), Error::Parquet { .. }));
        }

        #[test]
        fn test_parquet_write_permission_denied() {
            // Try to write to root directory (should fail with permission denied)
            let m = GenotypeMatrix::new(1);
            let result = m.into_parquet_file("/root/test.parquet");
            assert!(result.is_err());
            assert!(matches!(result.unwrap_err(), Error::Io { .. }));
        }
    }
}
