use arrow::array::*;

use arrow::record_batch::RecordBatch;
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

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GenotypeRecords {
    data: Vec<GenotypeRecord>,
    sort_status: u8, // 0: unsorted, 1: sorted by position, 2: sorted_by_genome
}

impl GenotypeRecords {
    // Create `GenotypeRecords` from a vector of `GenotypeRecord` and an indicator
    // of sort status.
    pub fn new(records: Vec<GenotypeRecord>, sort_status: u8) -> Self {
        Self {
            data: records,
            sort_status,
        }
    }

    pub fn merge(&mut self, other: Self) {
        self.data.extend(other.data);
        self.sort_status = 0;
    }

    pub fn sort_by_position(&mut self) {
        match self.sort_status {
            0 | 2 => {
                self.data.par_sort_unstable_by_key(|x| {
                    (x.get_pos_genome(), x.get_position(), x.get_allele())
                });
                self.sort_status = 1
            }
            1 => {}
            _ => {
                panic!("GenotypeRecords was constructed on a invalid sort_status value")
            }
        }
    }

    pub fn sort_by_genome(&mut self) {
        match self.sort_status {
            0 | 1 => {
                self.data.par_sort_unstable_by_key(|x| {
                    (x.get_genome_pos(), x.get_position(), x.get_allele())
                });
                self.sort_status = 2;
            }
            2 => {}
            _ => {
                panic!("GenotypeRecords was constructed on a invalid sort_status value")
            }
        }
    }

    pub fn is_sorted_by_postion(&self) -> bool {
        match self.sort_status {
            1 => true,
            0 | 2 => false,
            _ => {
                panic!("GenotypeRecords was constructed on a invalid sort_status value")
            }
        }
    }
    pub fn is_sorted_by_genome(&self) -> bool {
        match self.sort_status {
            2 => true,
            0 | 1 => false,
            _ => {
                panic!("GenotypeRecords was constructed on a invalid sort_status value")
            }
        }
    }

    pub fn iter_genome_pair_genotypes<'a>(
        &'a self,
        genome1: u32,
        genome2: u32,
    ) -> impl Iterator<Item = (u32, Option<u8>, Option<u8>)> + 'a {
        assert_eq!(self.sort_status, 2);
        let s1 = self.data.partition_point(|x| x.get_genome() < genome1);
        let e1 = self.data.partition_point(|x| x.get_genome() <= genome1);
        let s2 = self.data.partition_point(|x| x.get_genome() < genome2);
        let e2 = self.data.partition_point(|x| x.get_genome() <= genome2);
        let mergejoinby = self.data[s1..e1]
            .iter()
            .merge_join_by(self.data[s2..e2].iter(), |a, b| {
                (*a).get_position()
                    .partial_cmp(&(*b).get_position())
                    .unwrap()
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

    pub fn filter_multi_allelic_site(&mut self) {
        self.sort_by_position();
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
    }

    pub fn into_parquet_file(self, p: impl AsRef<Path>) {
        let v = self
            .data
            .into_iter()
            .map(|x| x.data.into_inner())
            .collect::<Vec<_>>();
        // array
        let u64values = UInt64Array::from(v);

        // save sort_status as field name
        let fieldname = format!("{}", self.sort_status);

        // record batch
        let batch =
            RecordBatch::try_from_iter(vec![(fieldname, Arc::new(u64values) as ArrayRef)]).unwrap();
        // writer
        let file = File::create(p.as_ref()).unwrap();
        // -- default writer properties
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();
        // write batch
        writer.write(&batch).expect("Writing batch");
        // writer must be closed to write footer
        writer.close().unwrap();
    }

    pub fn from_parquet_file(p: impl AsRef<Path>) -> Self {
        let file = File::open(p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        // get sort_status from the field name
        let sort_status = builder.schema().field(0).name().parse().unwrap();
        let mut reader = builder.build().unwrap();
        let mut records = Vec::<GenotypeRecord>::new();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            let rec_iter = record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt64Array>()
                .unwrap()
                .values()
                .iter()
                .map(|x| GenotypeRecord::new(*x));
            records.extend(rec_iter);
        }
        GenotypeRecords {
            data: records,
            sort_status,
        }
    }

    pub fn subset_by_genomes(&self, sorted_genome_ids: &[u32]) -> Result<Self, &'static str> {
        if !sorted_genome_ids
            .iter()
            .zip(sorted_genome_ids.iter().skip(1))
            .all(|(a, b)| *a <= *b)
        {
            return Err("the genome ids in `sorted_genome_ids` are not sorted ");
        }
        if !self.is_sorted_by_genome() {
            return Err("genotype are not sorted by genomes");
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
