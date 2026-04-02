use crate::error::{IshareError, Result};
use error_stack::*;

use arrow_array::{ArrayRef, RecordBatch, UInt64Array};

use bitvec::prelude::*;
use itertools::{
    EitherOrBoth::{self, Both, Left, Right},
    Itertools,
};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use slice_group_by::GroupBy;
use std::fmt::Debug;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use crate::traits::TotalOrd;

#[derive(Serialize, Deserialize, Clone)]
pub struct GenotypeRecord {
    data: BitArray<u64, Lsb0>,
    // pos: 32 bits, max value: 4,294,967,296 - 1 (4.2 billion)
    // (haploid) genome/sample id: 24 bits, max value: 16,777,216 - 1  (16.7 million)
    // Allele id: 8 bits, max value: 256 - 1
}
impl Debug for GenotypeRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "position={}, genome={}, allele={}",
            self.get_position(),
            self.get_genome(),
            self.get_allele()
        )
    }
}

impl GenotypeRecord {
    pub fn new(value: u64) -> Self {
        GenotypeRecord {
            data: BitArray::<u64, Lsb0>::new(value),
        }
    }
    pub fn set_sentinel(&mut self) {
        self.data.data = u64::MAX;
    }
    pub fn is_sentinel(&self) -> bool {
        self.data.data == u64::MAX
    }
    pub fn set_position(&mut self, value: u32) {
        // assert!((value >> 32) == 0);
        self.data[0..32].store(value);
    }
    pub fn get_position(&self) -> u32 {
        self.data[0..32].load()
    }
    pub fn set_genome(&mut self, value: u32) {
        assert!((value >> 24) == 0);
        self.data[32..56].store(value);
    }
    pub fn get_genome(&self) -> u32 {
        self.data[32..56].load()
    }
    pub fn set_allele(&mut self, value: u8) {
        //assert!((value >> 8) == 0);
        self.data[56..64].store(value);
    }
    pub fn get_allele(&self) -> u8 {
        self.data[56..64].load()
    }
    pub fn set(&mut self, position: u32, genome: u32, allele: u8) {
        self.set_position(position);
        self.set_genome(genome);
        self.set_allele(allele);
    }
    pub fn get(&self) -> (u32, u32, u8) {
        (self.get_position(), self.get_genome(), self.get_allele())
    }
    pub fn get_genome_pos(&self) -> (u32, u32) {
        (self.get_genome(), self.get_position())
    }
    pub fn get_pos_genome(&self) -> (u32, u32) {
        (self.get_position(), self.get_genome())
    }
    pub fn get_pos_allele(&self) -> (u32, u8) {
        (self.get_position(), self.get_allele())
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeRecordSortStatus {
    Unsorted,
    SortedByPositionGenomeAllele,
    SortedByGenomePositionAllele,
    SortedByPositionAlleleGenome,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GenotypeRecords {
    data: Vec<GenotypeRecord>,
    sort_status: GenotypeRecordSortStatus,
}

impl GenotypeRecords {
    // Create `GenotypeRecords` from a vector of `GenotypeRecord` and an indicator
    // of sort status.
    pub fn new(records: Vec<GenotypeRecord>, sort_status: GenotypeRecordSortStatus) -> Self {
        Self {
            data: records,
            sort_status,
        }
    }

    pub fn merge(&mut self, other: Self) {
        self.data.extend(other.data);
        self.sort_status = GenotypeRecordSortStatus::Unsorted;
    }

    pub fn assume_unsorted(&mut self) {
        self.sort_status = GenotypeRecordSortStatus::Unsorted;
    }

    /// useful for write genotype to vcf or genotype matrix
    pub fn sort_by_position_genome_allele(&mut self) -> Result<()> {
        match self.sort_status {
            GenotypeRecordSortStatus::SortedByGenomePositionAllele => Ok(()),
            _ => {
                self.data
                    .par_sort_unstable_by_key(|x| (x.get_pos_genome(), x.get_allele()));
                self.sort_status = GenotypeRecordSortStatus::SortedByPositionGenomeAllele;
                Ok(())
            }
        }
    }

    /// useful for counting shared rare variants by genome paris along chromosomal positions
    pub fn sort_by_genome_position_allele(&mut self) -> Result<()> {
        match self.sort_status {
            GenotypeRecordSortStatus::SortedByGenomePositionAllele => Ok(()),
            _ => {
                self.data
                    .par_sort_unstable_by_key(|x| (x.get_genome_pos(), x.get_allele()));
                self.sort_status = GenotypeRecordSortStatus::SortedByGenomePositionAllele;
                Ok(())
            }
        }
    }

    /// useful for count allele freuquency
    pub fn sort_by_position_allele_genome(&mut self) -> Result<()> {
        match self.sort_status {
            GenotypeRecordSortStatus::SortedByPositionAlleleGenome => Ok(()),
            _ => {
                self.data
                    .par_sort_unstable_by_key(|x| (x.get_pos_allele(), x.get_genome()));
                self.sort_status = GenotypeRecordSortStatus::SortedByPositionAlleleGenome;
                Ok(())
            }
        }
    }
    pub fn is_sorted_by_postion_genome_allele(&self) -> bool {
        matches!(
            self.sort_status,
            GenotypeRecordSortStatus::SortedByPositionGenomeAllele
        )
    }

    pub fn is_sorted_by_genome_position_allele_slow(&self) -> bool {
        let recs = self.records();
        recs.iter().zip(recs.iter().skip(1)).all(|(a, b)| {
            (a.get_genome_pos(), a.get_allele()) < (b.get_genome_pos(), b.get_allele())
        })
    }

    pub fn is_sorted_by_position_genome_allele_slow(&self) -> bool {
        let recs = self.records();
        recs.iter().zip(recs.iter().skip(1)).all(|(a, b)| {
            (a.get_pos_genome(), a.get_allele()) < (b.get_pos_genome(), b.get_allele())
        })
    }

    pub fn is_sorted_by_position_allele_genome_slow(&self) -> bool {
        let recs = self.records();
        recs.iter().zip(recs.iter().skip(1)).all(|(a, b)| {
            (a.get_pos_allele(), a.get_genome()) < (b.get_pos_allele(), b.get_genome())
        })
    }

    pub fn is_unsorted(&self) -> bool {
        matches!(self.sort_status, GenotypeRecordSortStatus::Unsorted)
    }

    pub fn is_sorted_by_postion_allele_genome(&self) -> bool {
        matches!(
            self.sort_status,
            GenotypeRecordSortStatus::SortedByPositionAlleleGenome
        )
    }

    pub fn get_sort_status(&self) -> GenotypeRecordSortStatus {
        self.sort_status
    }
    pub fn is_sorted_by_genome_position_allele(&self) -> bool {
        matches!(
            self.sort_status,
            GenotypeRecordSortStatus::SortedByGenomePositionAllele
        )
    }

    /// assume data is sorted by genome
    pub fn iterate_rv_shared_for_genome_pair(
        &self,
        genome1: u32,
        genome2: u32,
    ) -> impl Iterator<Item = (u32, u8)> + '_ {
        let s1 = self.data.partition_point(|x| x.get_genome() < genome1);
        let e1 = self.data.partition_point(|x| x.get_genome() <= genome1);
        let s2 = self.data.partition_point(|x| x.get_genome() < genome2);
        let e2 = self.data.partition_point(|x| x.get_genome() <= genome2);
        self.data[s1..e1]
            .iter()
            .merge_join_by(self.data[s2..e2].iter(), |a, b| {
                (*a).get_position().total_cmp(&(*b).get_position())
            })
            .filter_map(move |res| match res {
                itertools::EitherOrBoth::Both(a, b) => {
                    if a.get_allele() == b.get_allele() {
                        Some((a.get_position(), a.get_allele()))
                    } else {
                        None
                    }
                }
                _ => None,
            })
    }

    /// assert data is sorted by position and allele
    pub fn get_allele_count_by_position_allele(&self, pos: u32, allele: u8) -> usize {
        let s = self
            .records()
            .partition_point(|r| r.get_pos_allele() < (pos, allele));
        let e = self
            .records()
            .partition_point(|r| r.get_pos_allele() <= (pos, allele));
        e - s
    }

    pub fn iter_genome_pair_genotypes(
        &self,
        genome1: u32,
        genome2: u32,
    ) -> impl Iterator<Item = (u32, Option<u8>, Option<u8>)> + '_ {
        assert!(self.is_sorted_by_genome_position_allele());
        let s1 = self.data.partition_point(|x| x.get_genome() < genome1);
        let e1 = self.data.partition_point(|x| x.get_genome() <= genome1);
        let s2 = self.data.partition_point(|x| x.get_genome() < genome2);
        let e2 = self.data.partition_point(|x| x.get_genome() <= genome2);
        let mergejoinby = self.data[s1..e1]
            .iter()
            .merge_join_by(self.data[s2..e2].iter(), |a, b| {
                (*a).get_position().total_cmp(&(*b).get_position())
            })
            .map(|res| match res {
                Both(a, b) => (a.get_position(), Some(a.get_allele()), Some(b.get_allele())),
                Left(a) => (a.get_position(), Some(a.get_allele()), None),
                Right(b) => (b.get_position(), None, Some(b.get_allele())),
            });
        mergejoinby
    }

    pub fn records(&self) -> &[GenotypeRecord] {
        &self.data
    }
    pub fn records_mut(&mut self) -> &mut Vec<GenotypeRecord> {
        &mut self.data
    }

    pub fn filter_multi_allelic_site(&mut self) -> Result<()> {
        self.sort_by_position_genome_allele()?;
        use slice_group_by::*;
        for blk in self
            .data
            .as_mut_slice()
            .linear_group_by_key_mut(|r| r.get_position())
        {
            // identify position with multile alleles
            let allele = blk[0].get_allele();
            if blk.iter().any(|x| x.get_allele() != allele) {
                // mark for all records at this position for deletion
                blk.iter_mut().for_each(|r| r.set_position(u32::MAX));
            }
        }
        self.data.retain(|x| x.get_position() != u32::MAX);
        Ok(())
    }

    pub fn into_parquet_file(self, p: impl AsRef<Path>) -> Result<()> {
        let v = self
            .data
            .into_iter()
            .map(|x| x.data.into_inner())
            .collect::<Vec<_>>();
        // array
        let u64values = UInt64Array::from(v);

        // save sort_status as field name
        let fieldname = format!("{:?}", self.sort_status);

        // record batch
        let batch = RecordBatch::try_from_iter(vec![(fieldname, Arc::new(u64values) as ArrayRef)])
            .change_context(IshareError::RareGenotype)?;
        // writer
        let file = File::create(p.as_ref()).change_context(IshareError::RareGenotype)?;
        // -- default writer properties
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))
            .change_context(IshareError::RareGenotype)?;
        // write batch
        writer
            .write(&batch)
            .change_context(IshareError::RareGenotype)?;
        // writer must be closed to write footer
        writer.close().change_context(IshareError::RareGenotype)?;
        Ok(())
    }

    pub fn from_parquet_file(p: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(p).change_context(IshareError::RareGenotype)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .change_context(IshareError::RareGenotype)?;
        // get sort_status from the field name
        let sort_status = match builder.schema().field(0).name().as_str() {
            "SortedByPositionGenomeAllele" => {
                GenotypeRecordSortStatus::SortedByPositionGenomeAllele
            }
            "SortedByGenomePositionAllele" => {
                GenotypeRecordSortStatus::SortedByGenomePositionAllele
            }
            "SortedByPositionAlleleGenome" => {
                GenotypeRecordSortStatus::SortedByPositionAlleleGenome
            }
            _ => GenotypeRecordSortStatus::Unsorted,
        };
        let mut reader = builder.build().change_context(IshareError::RareGenotype)?;
        let mut records = Vec::<GenotypeRecord>::new();
        for record_batch in &mut reader {
            let record_batch = record_batch.change_context(IshareError::RareGenotype)?;
            let rec_iter = record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .ok_or(IshareError::EmptyOption)?
                .values()
                .iter()
                .map(|x| GenotypeRecord::new(*x));
            records.extend(rec_iter);
        }
        Ok(GenotypeRecords {
            data: records,
            sort_status,
        })
    }

    pub fn subset_by_genomes(&self, sorted_genome_ids: &[u32]) -> Result<Self> {
        if !sorted_genome_ids
            .iter()
            .zip(sorted_genome_ids.iter().skip(1))
            .all(|(a, b)| *a <= *b)
        {
            bail!(IshareError::RuntimeCheck
                .into_report()
                .attach("GenomeIdsNotSorted"));
        }
        if !self.is_sorted_by_genome_position_allele() {
            bail!(IshareError::RuntimeCheck
                .into_report()
                .attach("GenomeIdsNotSorted"));
        }

        let mut v = vec![];
        for e in self
            .data
            .linear_group_by_key(|r| r.get_genome())
            .merge_join_by(sorted_genome_ids.linear_group(), |a, b| {
                a[0].get_genome().cmp(&b[0])
            })
        {
            match e {
                EitherOrBoth::Left(_) => {}
                EitherOrBoth::Right(_) => {}
                EitherOrBoth::Both(left, _right) => {
                    v.extend_from_slice(left);
                }
            }
        }
        Ok(GenotypeRecords {
            data: v,
            sort_status: self.sort_status,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    // GenotypeRecord Tests
    mod genotype_record_tests {
        use super::*;

        #[test]
        fn test_new_and_basic_operations() {
            let record = GenotypeRecord::new(0);

            // Test initial state
            assert_eq!(record.get_position(), 0);
            assert_eq!(record.get_genome(), 0);
            assert_eq!(record.get_allele(), 0);
            assert!(!record.is_sentinel());
        }

        #[test]
        fn test_bit_field_boundaries() {
            let mut record = GenotypeRecord::new(0);

            // Test position field (32 bits)
            record.set_position(u32::MAX);
            assert_eq!(record.get_position(), u32::MAX);

            // Test genome field (24 bits)
            let max_genome = (1u32 << 24) - 1; // 16,777,215
            record.set_genome(max_genome);
            assert_eq!(record.get_genome(), max_genome);

            // Test allele field (8 bits)
            record.set_allele(u8::MAX);
            assert_eq!(record.get_allele(), u8::MAX);
        }

        #[test]
        #[should_panic]
        fn test_genome_field_overflow() {
            let mut record = GenotypeRecord::new(0);
            // This should panic due to assertion
            let invalid_genome = 1u32 << 24; // 16,777,216 (too large for 24 bits)
            record.set_genome(invalid_genome);
        }

        #[test]
        fn test_set_and_get_operations() {
            let mut record = GenotypeRecord::new(0);

            // Test set method
            record.set(123456, 789, 42);
            let (pos, genome, allele) = record.get();
            assert_eq!(pos, 123456);
            assert_eq!(genome, 789);
            assert_eq!(allele, 42);
        }

        #[test]
        fn test_sentinel_operations() {
            let mut record = GenotypeRecord::new(0);

            assert!(!record.is_sentinel());
            record.set_sentinel();
            assert!(record.is_sentinel());

            // After setting sentinel, all bits should be 1
            assert_eq!(record.data.data, u64::MAX);
        }

        #[test]
        fn test_getter_combinations() {
            let mut record = GenotypeRecord::new(0);
            record.set(100, 200, 50);

            assert_eq!(record.get_genome_pos(), (200, 100));
            assert_eq!(record.get_pos_genome(), (100, 200));
            assert_eq!(record.get_pos_allele(), (100, 50));
        }

        #[test]
        fn test_bit_manipulation_precision() {
            let mut record = GenotypeRecord::new(0);

            // Set specific bit patterns to verify no bit bleeding
            record.set_position(0xAAAAAAAA);
            record.set_genome(0x555555); // 24-bit value
            record.set_allele(0xCC);

            assert_eq!(record.get_position(), 0xAAAAAAAA);
            assert_eq!(record.get_genome(), 0x555555);
            assert_eq!(record.get_allele(), 0xCC);
        }

        #[test]
        fn test_zero_values() {
            let mut record = GenotypeRecord::new(u64::MAX);
            record.set(0, 0, 0);

            assert_eq!(record.get_position(), 0);
            assert_eq!(record.get_genome(), 0);
            assert_eq!(record.get_allele(), 0);
        }
    }

    // GenotypeRecords Tests
    mod genotype_records_tests {
        use super::*;

        fn create_test_records() -> Vec<GenotypeRecord> {
            let mut records = Vec::new();

            // Create records with different positions and genomes for testing
            let test_data = [
                (100, 1, 1),
                (200, 1, 2),
                (100, 2, 1),
                (300, 2, 1),
                (150, 1, 1),
            ];

            for (pos, genome, allele) in test_data {
                let mut record = GenotypeRecord::new(0);
                record.set(pos, genome, allele);
                records.push(record);
            }

            records
        }

        #[test]
        fn test_new_and_basic_properties() {
            let records = create_test_records();
            let genotype_records =
                GenotypeRecords::new(records.clone(), GenotypeRecordSortStatus::Unsorted);

            assert_eq!(genotype_records.records().len(), records.len());
            assert!(genotype_records.is_unsorted());
        }

        #[test]
        fn test_merge_operations() {
            let records1 = create_test_records();
            let records2 = create_test_records();

            let mut genotype_records1 =
                GenotypeRecords::new(records1.clone(), GenotypeRecordSortStatus::Unsorted);
            let genotype_records2 =
                GenotypeRecords::new(records2.clone(), GenotypeRecordSortStatus::Unsorted);

            let original_len = genotype_records1.records().len();
            genotype_records1.merge(genotype_records2);

            assert_eq!(genotype_records1.records().len(), original_len * 2);
            assert!(genotype_records1.is_unsorted()); // Should reset to unsorted
        }

        #[test]
        fn test_sort_by_position() -> Result<()> {
            use GenotypeRecordSortStatus::*;
            let records = create_test_records();
            let mut genotype_records = GenotypeRecords::new(records, Unsorted);

            // Initial state should be unsorted
            assert!(!genotype_records.is_sorted_by_postion_genome_allele());

            // Sort by position
            genotype_records.sort_by_position_genome_allele()?;
            assert!(genotype_records.is_sorted_by_postion_genome_allele());
            assert!(!genotype_records.is_sorted_by_genome_position_allele());

            // Verify sorting correctness
            let positions: Vec<u32> = genotype_records
                .records()
                .iter()
                .map(|r| r.get_position())
                .collect();
            let mut sorted_positions = positions.clone();
            sorted_positions.sort();
            assert_eq!(positions, sorted_positions);

            Ok(())
        }

        #[test]
        fn test_sort_by_genome() -> Result<()> {
            use GenotypeRecordSortStatus::*;
            let records = create_test_records();
            let mut genotype_records = GenotypeRecords::new(records, Unsorted);

            // Sort by genome
            genotype_records.sort_by_genome_position_allele()?;
            assert!(genotype_records.is_sorted_by_genome_position_allele());
            assert!(!genotype_records.is_sorted_by_postion_genome_allele());

            // Verify sorting correctness
            let genomes: Vec<u32> = genotype_records
                .records()
                .iter()
                .map(|r| r.get_genome())
                .collect();
            let mut sorted_genomes = genomes.clone();
            sorted_genomes.sort();
            assert_eq!(genomes, sorted_genomes);

            Ok(())
        }

        #[test]
        fn test_sort_status_transitions() -> Result<()> {
            use GenotypeRecordSortStatus::*;
            let records = create_test_records();
            let mut genotype_records = GenotypeRecords::new(records, Unsorted);

            // Unsorted -> Position sorted
            genotype_records.sort_by_position_genome_allele()?;
            assert!(genotype_records.is_sorted_by_postion_genome_allele());

            // Position sorted -> Genome sorted
            genotype_records.sort_by_genome_position_allele()?;
            assert!(genotype_records.is_sorted_by_genome_position_allele());

            // Genome sorted -> Position sorted
            genotype_records.sort_by_position_genome_allele()?;
            assert!(genotype_records.is_sorted_by_postion_genome_allele());

            Ok(())
        }

        #[test]
        fn test_iter_genome_pair_genotypes() -> Result<()> {
            // Create specific test data for genome pair iteration
            let mut records = Vec::new();

            // Genome 1 has variants at positions 100, 200, 400
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(100, 1, 1);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(200, 1, 2);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(400, 1, 1);
                r
            });

            // Genome 2 has variants at positions 150, 200, 300
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(150, 2, 1);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(200, 2, 3);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(300, 2, 2);
                r
            });

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            genotype_records.sort_by_genome_position_allele()?;

            let pairs: Vec<_> = genotype_records.iter_genome_pair_genotypes(1, 2).collect();

            // Expected: (position, genome1_allele, genome2_allele)
            assert_eq!(
                pairs,
                vec![
                    (100, Some(1), None),    // Only genome 1
                    (150, None, Some(1)),    // Only genome 2
                    (200, Some(2), Some(3)), // Both genomes
                    (300, None, Some(2)),    // Only genome 2
                    (400, Some(1), None),    // Only genome 1
                ]
            );

            Ok(())
        }

        #[test]
        fn test_filter_multi_allelic_site() -> Result<()> {
            let mut records = Vec::new();

            // Position 100: single allele (should keep)
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(100, 1, 1);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(100, 2, 1);
                r
            });

            // Position 200: multiple alleles (should remove)
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(200, 1, 1);
                r
            });
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(200, 2, 2);
                r
            }); // Different allele

            // Position 300: single allele (should keep)
            records.push({
                let mut r = GenotypeRecord::new(0);
                r.set(300, 1, 3);
                r
            });

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            let original_count = genotype_records.records().len();

            genotype_records.filter_multi_allelic_site()?;

            // Should keep positions 100 and 300, remove position 200
            assert_eq!(genotype_records.records().len(), 3);
            assert!(original_count > genotype_records.records().len());

            // Verify no records at position 200 remain
            assert!(!genotype_records
                .records()
                .iter()
                .any(|r| r.get_position() == 200));

            Ok(())
        }

        #[test]
        fn test_subset_by_genomes() -> Result<()> {
            let mut records = Vec::new();

            // Create records for genomes 1, 2, 3, 5
            for genome in [1, 2, 3, 5] {
                for pos in [100, 200] {
                    let mut record = GenotypeRecord::new(0);
                    record.set(pos, genome, 1);
                    records.push(record);
                }
            }

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            genotype_records.sort_by_genome_position_allele()?;

            // Subset to genomes [2, 5] (must be sorted)
            let subset_genomes = vec![2, 5];
            let subset = genotype_records.subset_by_genomes(&subset_genomes)?;

            assert_eq!(subset.records().len(), 4); // 2 genomes × 2 positions

            // Verify all records are from requested genomes
            for record in subset.records() {
                assert!(subset_genomes.contains(&record.get_genome()));
            }

            Ok(())
        }

        #[test]
        fn test_subset_by_genomes_requires_sorted_input() {
            let records = create_test_records();
            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            genotype_records.sort_by_genome_position_allele().unwrap();

            // Unsorted genome list should fail
            let unsorted_genomes = vec![3, 1, 2];
            assert!(genotype_records
                .subset_by_genomes(&unsorted_genomes)
                .is_err());
        }

        #[test]
        fn test_subset_by_genomes_requires_genome_sorted_records() {
            let records = create_test_records();
            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted); // Position sorted, not genome sorted

            genotype_records.sort_by_position_genome_allele().unwrap();

            let sorted_genomes = vec![1, 2];
            assert!(genotype_records.subset_by_genomes(&sorted_genomes).is_err());
        }

        #[test]
        fn test_empty_records() -> Result<()> {
            let empty_records = GenotypeRecords::new(vec![], GenotypeRecordSortStatus::Unsorted);

            assert_eq!(empty_records.records().len(), 0);
            assert!(!empty_records.is_sorted_by_postion_genome_allele());
            assert!(!empty_records.is_sorted_by_genome_position_allele());

            // Operations on empty records should work
            let pairs: Vec<_> = {
                let mut temp = empty_records.clone();
                temp.sort_by_genome_position_allele()?;
                temp.iter_genome_pair_genotypes(1, 2).collect()
            };
            assert_eq!(pairs.len(), 0);

            Ok(())
        }
    }

    // File I/O and Serialization Tests
    mod serialization_tests {
        use super::*;

        fn create_test_data() -> GenotypeRecords {
            let mut records = Vec::new();

            // Create diverse test data
            let test_cases = [
                (1000, 1, 1),
                (2000, 1, 2),
                (1500, 2, 1),
                (3000, 3, 255),                    // Max allele value
                (u32::MAX, (1u32 << 24) - 1, 128), // Max position and genome values
            ];

            for (pos, genome, allele) in test_cases {
                let mut record = GenotypeRecord::new(0);
                record.set(pos, genome, allele);
                records.push(record);
            }

            GenotypeRecords::new(
                records,
                GenotypeRecordSortStatus::SortedByPositionGenomeAllele,
            ) // Position sorted
        }

        #[test]
        fn test_parquet_round_trip() -> Result<()> {
            let original_data = create_test_data();
            let temp_dir = tempdir().unwrap();
            let temp_path = temp_dir.path().join("test.parquet");

            // Write to parquet
            original_data.clone().into_parquet_file(&temp_path)?;

            // Read back from parquet
            let loaded_data = GenotypeRecords::from_parquet_file(&temp_path)?;

            // Verify data integrity
            assert_eq!(loaded_data.sort_status, original_data.sort_status);
            assert_eq!(loaded_data.records().len(), original_data.records().len());

            for (original, loaded) in original_data
                .records()
                .iter()
                .zip(loaded_data.records().iter())
            {
                assert_eq!(original.get(), loaded.get());
            }

            Ok(())
        }

        #[test]
        fn test_parquet_preserves_sort_status() -> Result<()> {
            for sort_status in [
                GenotypeRecordSortStatus::Unsorted,
                GenotypeRecordSortStatus::SortedByPositionGenomeAllele,
                GenotypeRecordSortStatus::SortedByGenomePositionAllele,
                GenotypeRecordSortStatus::SortedByPositionAlleleGenome,
            ] {
                let mut records = create_test_data();
                records.sort_status = sort_status;

                let temp_dir = tempdir().unwrap();
                let temp_path = temp_dir
                    .path()
                    .join(format!("test_{sort_status:?}.parquet"));

                records.clone().into_parquet_file(&temp_path)?;
                let loaded = GenotypeRecords::from_parquet_file(&temp_path)?;

                assert_eq!(loaded.sort_status, sort_status);
            }

            Ok(())
        }

        #[test]
        fn test_parquet_empty_data() -> Result<()> {
            let empty_data = GenotypeRecords::new(
                vec![],
                GenotypeRecordSortStatus::SortedByGenomePositionAllele,
            );
            let temp_dir = tempdir().unwrap();
            let temp_path = temp_dir.path().join("empty.parquet");

            empty_data.clone().into_parquet_file(&temp_path)?;
            let loaded = GenotypeRecords::from_parquet_file(&temp_path)?;

            assert_eq!(loaded.records().len(), 0);
            assert_eq!(
                loaded.sort_status,
                GenotypeRecordSortStatus::SortedByGenomePositionAllele
            );

            Ok(())
        }

        #[test]
        fn test_parquet_large_dataset() -> Result<()> {
            // Create larger dataset to test performance and correctness
            let mut records = Vec::new();
            for i in 0..10000 {
                let mut record = GenotypeRecord::new(0);
                record.set(i * 100, i % 1000, (i % 256) as u8);
                records.push(record);
            }

            let original_data = GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            let temp_dir = tempdir().unwrap();
            let temp_path = temp_dir.path().join("large.parquet");

            original_data.clone().into_parquet_file(&temp_path)?;
            let loaded_data = GenotypeRecords::from_parquet_file(&temp_path)?;

            assert_eq!(loaded_data.records().len(), original_data.records().len());
            assert_eq!(loaded_data.sort_status, original_data.sort_status);

            Ok(())
        }

        #[test]
        fn test_parquet_file_not_found() {
            let result = GenotypeRecords::from_parquet_file("/nonexistent/path.parquet");
            assert!(result.is_err());
        }

        #[test]
        fn test_parquet_invalid_directory() {
            let invalid_path = "/invalid/directory/test.parquet";
            let data = create_test_data();
            let result = data.into_parquet_file(invalid_path);
            assert!(result.is_err());
        }
    }

    // Edge Cases and Error Handling Tests
    mod edge_cases_tests {
        use super::*;

        #[test]
        fn test_extreme_values() {
            let mut record = GenotypeRecord::new(0);

            // Test maximum values for each field
            record.set_position(u32::MAX);
            record.set_genome((1u32 << 24) - 1); // 24-bit max
            record.set_allele(u8::MAX);

            assert_eq!(record.get_position(), u32::MAX);
            assert_eq!(record.get_genome(), (1u32 << 24) - 1);
            assert_eq!(record.get_allele(), u8::MAX);
        }

        #[test]
        fn test_bit_field_isolation() {
            let mut record = GenotypeRecord::new(0);

            // Set each field to maximum, verify others unaffected
            record.set_position(u32::MAX);
            assert_eq!(record.get_genome(), 0);
            assert_eq!(record.get_allele(), 0);

            record.set_genome((1u32 << 24) - 1);
            assert_eq!(record.get_position(), u32::MAX);
            assert_eq!(record.get_allele(), 0);

            record.set_allele(u8::MAX);
            assert_eq!(record.get_position(), u32::MAX);
            assert_eq!(record.get_genome(), (1u32 << 24) - 1);
        }

        #[test]
        fn test_single_record_operations() -> Result<()> {
            let record = GenotypeRecord::new(0);
            let mut record = record;
            record.set(12345, 678, 90);

            let mut records =
                GenotypeRecords::new(vec![record], GenotypeRecordSortStatus::Unsorted);

            // All operations should work with single record
            records.sort_by_position_genome_allele()?;
            assert!(records.is_sorted_by_postion_genome_allele());

            records.sort_by_genome_position_allele()?;
            assert!(records.is_sorted_by_genome_position_allele());

            records.filter_multi_allelic_site()?;
            assert_eq!(records.records().len(), 1);

            // Need to sort by genome again since filter_multi_allelic_site calls sort_by_position internally
            records.sort_by_genome_position_allele()?;
            let subset = records.subset_by_genomes(&[678])?;
            assert_eq!(subset.records().len(), 1);

            Ok(())
        }

        #[test]
        fn test_duplicate_positions_sorting() -> Result<()> {
            let mut records = Vec::new();

            // Create multiple records at same position
            for genome in 1..=5 {
                let mut record = GenotypeRecord::new(0);
                record.set(100, genome, 1);
                records.push(record);
            }

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            genotype_records.sort_by_position_genome_allele()?;

            // All should be at position 100
            for record in genotype_records.records() {
                assert_eq!(record.get_position(), 100);
            }

            Ok(())
        }

        #[test]
        fn test_memory_efficiency_patterns() {
            // Test that bit packing is working efficiently
            let mut record1 = GenotypeRecord::new(0);
            let mut record2 = GenotypeRecord::new(0);

            record1.set(1, 1, 1);
            record2.set(2, 2, 2);

            // Records should have different underlying bit patterns
            assert_ne!(record1.data.data, record2.data.data);

            // But similar records should be identical
            let mut record3 = GenotypeRecord::new(0);
            record3.set(1, 1, 1);
            assert_eq!(record1.data.data, record3.data.data);
        }

        #[test]
        fn test_large_genome_ids() -> Result<()> {
            let mut records = Vec::new();

            // Test with large valid genome IDs (close to 24-bit limit)
            let large_genome_ids = [
                (1u32 << 23) - 1, // Half of 24-bit space
                (1u32 << 23),     // Middle of 24-bit space
                (1u32 << 24) - 2, // Near maximum
                (1u32 << 24) - 1, // Maximum valid value
            ];

            for (i, genome_id) in large_genome_ids.iter().enumerate() {
                let mut record = GenotypeRecord::new(0);
                record.set(i as u32 * 1000, *genome_id, (i + 1) as u8);
                records.push(record);
            }

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);
            genotype_records.sort_by_genome_position_allele()?;

            // Verify all genome IDs preserved correctly
            for (record, expected_genome) in genotype_records
                .records()
                .iter()
                .zip(large_genome_ids.iter())
            {
                assert_eq!(record.get_genome(), *expected_genome);
            }

            Ok(())
        }
    }

    // Performance and Parallel Processing Tests
    mod performance_tests {
        use super::*;

        #[test]
        fn test_parallel_sort_consistency() -> Result<()> {
            // Create large dataset to trigger parallel sorting
            let mut records = Vec::new();
            use std::collections::hash_map::DefaultHasher;
            use std::hash::{Hash, Hasher};

            // Generate pseudo-random but deterministic test data
            for i in 0..50000 {
                let mut hasher = DefaultHasher::new();
                i.hash(&mut hasher);
                let hash = hasher.finish();

                let mut record = GenotypeRecord::new(0);
                record.set(
                    (hash % (u32::MAX as u64)) as u32,
                    ((hash >> 32) % ((1u64 << 24) - 1)) as u32,
                    ((hash >> 40) % 256) as u8,
                );
                records.push(record);
            }

            let mut genotype_records1 =
                GenotypeRecords::new(records.clone(), GenotypeRecordSortStatus::Unsorted);
            let mut genotype_records2 =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);

            // Both should produce identical results despite parallel processing
            genotype_records1.sort_by_position_genome_allele()?;
            genotype_records2.sort_by_position_genome_allele()?;

            assert_eq!(
                genotype_records1.records().len(),
                genotype_records2.records().len()
            );
            for (r1, r2) in genotype_records1
                .records()
                .iter()
                .zip(genotype_records2.records().iter())
            {
                assert_eq!(r1.get(), r2.get());
            }

            Ok(())
        }

        #[test]
        fn test_large_dataset_operations() -> Result<()> {
            // Test with moderately large dataset
            let mut records = Vec::new();
            for i in 0..25000 {
                let mut record = GenotypeRecord::new(0);
                record.set(i * 2, i % 1000, ((i % 255) + 1) as u8);
                records.push(record);
            }

            let mut genotype_records =
                GenotypeRecords::new(records, GenotypeRecordSortStatus::Unsorted);

            // Test all major operations work with large datasets
            genotype_records.sort_by_position_genome_allele()?;
            assert!(genotype_records.is_sorted_by_postion_genome_allele());

            genotype_records.sort_by_genome_position_allele()?;
            assert!(genotype_records.is_sorted_by_genome_position_allele());

            // Filter multi-allelic sites needs position sorted data
            genotype_records.filter_multi_allelic_site()?;

            // After filtering, it should still be position sorted
            assert!(genotype_records.is_sorted_by_postion_genome_allele());

            // Verify data integrity maintained
            assert!(!genotype_records.records().is_empty());

            Ok(())
        }
    }
}
