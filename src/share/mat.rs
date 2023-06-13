use crate::io::{FromParquet, IntoParquet};
use arrow::array::*;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;
use std::{collections::HashMap, fmt::Debug};

#[derive(Debug, Clone)]
pub struct ResultMatrix {
    row_genomes: Vec<u32>,
    row_genomes_map: HashMap<u32, u32>,
    col_genomes: Vec<u32>,
    col_genomes_map: HashMap<u32, u32>,
    // parquet is very slow when using f32 type but very fast using integer
    // the value in the matrix is  (jaccard value) * 1e6 as u32
    data: Vec<u32>,
}

impl ResultMatrix {
    pub fn new_from_shape(nrow: u32, ncol: u32) -> Self {
        let row_genomes: Vec<u32> = (0..nrow).collect();
        let col_genomes: Vec<u32> = (0..ncol).collect();
        Self::new(row_genomes, col_genomes)
    }
    pub fn new(mut row_genomes: Vec<u32>, mut col_genomes: Vec<u32>) -> Self {
        row_genomes.sort();
        col_genomes.sort();

        // check they are unique
        assert!(row_genomes
            .iter()
            .zip(row_genomes.iter().skip(1))
            .all(|(a, b)| *a < *b));
        assert!(col_genomes
            .iter()
            .zip(col_genomes.iter().skip(1))
            .all(|(a, b)| *a < *b));

        // build hash map
        let row_genomes_map: HashMap<_, _> = row_genomes
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i as u32))
            .collect();
        let col_genomes_map: HashMap<_, _> = col_genomes
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i as u32))
            .collect();

        // allocation and initiation
        let n = row_genomes.len() * col_genomes.len();
        let mut data = Vec::new();
        data.resize(n, 0);

        Self {
            row_genomes,
            row_genomes_map,
            col_genomes,
            col_genomes_map,
            data,
        }
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.row_genomes.len(), self.col_genomes.len())
    }

    pub fn set_at(&mut self, row_idx: u32, col_idx: u32, v: u32) {
        let idx = (row_idx as usize) * self.col_genomes.len() + (col_idx as usize);
        self.data[idx] = v;
    }
    pub fn set_at_genomes(&mut self, row_genome: u32, col_genome: u32, v: u32) {
        let row_idx = self.row_genomes_map[&row_genome];
        let col_idx = self.col_genomes_map[&col_genome];
        let idx = (row_idx as usize) * self.col_genomes.len() + (col_idx as usize);
        self.data[idx] = v;
    }

    pub fn get_at(&self, row_idx: u32, col_idx: u32) -> u32 {
        let idx = (row_idx as usize) * self.col_genomes.len() + (col_idx as usize);
        self.data[idx]
    }

    pub fn get_at_genomes(&mut self, row_genome: u32, col_genome: u32) -> u32 {
        let row_idx = self.row_genomes_map[&row_genome];
        let col_idx = self.col_genomes_map[&col_genome];
        let idx = (row_idx as usize) * self.col_genomes.len() + (col_idx as usize);
        self.data[idx]
    }

    pub fn contains_matrix_by_genomes(&self, other: &Self) -> bool {
        let row_is_subset = other
            .row_genomes
            .iter()
            .all(|x| self.row_genomes_map.contains_key(x));
        let col_is_subset = other
            .col_genomes
            .iter()
            .all(|x| self.col_genomes_map.contains_key(x));
        row_is_subset && col_is_subset
    }

    pub fn update_from(&mut self, other: &Self) {
        assert!(self.contains_matrix_by_genomes(other));
        let (nrow, ncol) = other.shape();
        for i in 0..nrow {
            for j in 0..ncol {
                let v = other.get_at(i as u32, j as u32);

                let row_idx = self.row_genomes_map[&other.row_genomes[i]];
                let col_idx = self.col_genomes_map[&other.col_genomes[j]];

                self.set_at(row_idx, col_idx, v);
            }
        }
    }

    pub fn is_equal(&self, other: &Self) -> bool {
        (self.row_genomes == other.row_genomes)
            && (self.col_genomes == other.col_genomes)
            && (self.data == other.data)
    }
}

impl IntoParquet for ResultMatrix {
    fn into_parquet(&mut self, p: impl AsRef<Path>) {
        use std::mem::take;
        // build array from genotype matrix
        let mat = Arc::new(UInt32Array::from(take(&mut self.data))) as ArrayRef;
        let row = Arc::new(UInt32Array::from(take(&mut self.row_genomes))) as ArrayRef;
        let col = Arc::new(UInt32Array::from(take(&mut self.col_genomes))) as ArrayRef;

        let write_array = |fieldname, arr, p: &Path| {
            let batch =
                RecordBatch::try_from_iter(vec![(fieldname, Arc::new(arr) as ArrayRef)]).unwrap();
            // writer
            // use std::io::BufWriter;
            let file = BufWriter::new(File::create(p).unwrap());
            // -- default writer properties
            let props = WriterProperties::builder().build();
            let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();
            // write batch
            writer.write(&batch).expect("Writing batch");
            // writer must be closed to write footer
            writer.close().unwrap();
        };

        let row_file = p.as_ref().with_extension("row");
        write_array("row", row, &row_file);

        let col_file = p.as_ref().with_extension("col");
        write_array("col", col, &col_file);

        write_array("jaccard", mat, p.as_ref());
    }
}

impl FromParquet for ResultMatrix {
    fn from_parquet(p: impl AsRef<Path>) -> Self {
        let mut row_genomes = Vec::<u32>::new();
        let mut col_genomes = Vec::<u32>::new();
        let mut data = Vec::<u32>::new();

        let file = File::open(&p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .into_iter()
                .for_each(|x| data.push(x.unwrap()));
        }

        let file = File::open(p.as_ref().with_extension("row")).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .into_iter()
                .for_each(|x| row_genomes.push(x.unwrap()));
        }

        let file = File::open(p.as_ref().with_extension("col")).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap()
                .into_iter()
                .for_each(|x| col_genomes.push(x.unwrap()));
        }

        let mut matrix = Self::new(row_genomes, col_genomes);
        matrix.data = data;

        matrix
    }
}
