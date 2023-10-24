use crate::io::{FromArrowArray, FromParquet, IntoArrowArray, IntoParquet};
use ahash::HashMap;
use arrow::array::*;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::fmt::Debug;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct NamedMatrix<T>
where
    T: Copy + Default + Sized + PartialEq,
    Vec<T>: IntoArrowArray,
    [T]: FromArrowArray,
{
    row_names: Vec<u32>,
    row_names_map: HashMap<u32, u32>,
    col_names: Vec<u32>,
    col_names_map: HashMap<u32, u32>,
    // parquet is very slow when using f32 type but very fast using integer
    // the value in the matrix is  (jaccard value) * 1e6 as u32
    data: Vec<T>,
}

impl<T> NamedMatrix<T>
where
    T: Copy + Default + Sized + PartialEq,
    Vec<T>: IntoArrowArray,
    [T]: FromArrowArray,
{
    pub fn new_from_shape(nrow: u32, ncol: u32) -> Self {
        let row_names: Vec<u32> = (0..nrow).collect();
        let col_names: Vec<u32> = (0..ncol).collect();
        Self::new(row_names, col_names)
    }
    pub fn new_from_shape_and_data(nrow: u32, ncol: u32, data: Vec<T>) -> Self {
        let row_names: Vec<u32> = (0..nrow).collect();
        let col_names: Vec<u32> = (0..ncol).collect();
        let row_names_map: HashMap<u32, u32> = row_names.iter().map(|x| (*x, *x)).collect();
        let col_names_map: HashMap<u32, u32> = col_names.iter().map(|x| (*x, *x)).collect();
        assert_eq!(nrow as usize * ncol as usize, data.len());
        Self {
            row_names,
            row_names_map,
            col_names,
            col_names_map,
            data,
        }
    }
    pub fn new(mut row_names: Vec<u32>, mut col_names: Vec<u32>) -> Self {
        row_names.sort();
        col_names.sort();

        // check they are unique
        assert!(row_names
            .iter()
            .zip(row_names.iter().skip(1))
            .all(|(a, b)| *a < *b));
        assert!(col_names
            .iter()
            .zip(col_names.iter().skip(1))
            .all(|(a, b)| *a < *b));

        // build hash map
        let row_names_map: HashMap<_, _> = row_names
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i as u32))
            .collect();
        let col_names_map: HashMap<_, _> = col_names
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i as u32))
            .collect();

        // allocation and initiation
        let n = row_names.len() * col_names.len();
        let mut data = Vec::new();
        data.resize(n, T::default());

        Self {
            row_names,
            row_names_map,
            col_names,
            col_names_map,
            data,
        }
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.row_names.len(), self.col_names.len())
    }

    pub fn set_by_positions(&mut self, row_idx: u32, col_idx: u32, v: T) {
        let idx = (row_idx as usize) * self.col_names.len() + (col_idx as usize);
        self.data[idx] = v;
    }
    pub fn set_by_names(&mut self, row_genome: u32, col_genome: u32, v: T) {
        let row_idx = self.row_names_map[&row_genome];
        let col_idx = self.col_names_map[&col_genome];
        let idx = (row_idx as usize) * self.col_names.len() + (col_idx as usize);
        self.data[idx] = v;
    }

    pub fn get_by_positions(&self, row_idx: u32, col_idx: u32) -> T {
        let idx = (row_idx as usize) * self.col_names.len() + (col_idx as usize);
        self.data[idx]
    }

    pub fn get_by_names(&mut self, row_genome: u32, col_genome: u32) -> T {
        let row_idx = self.row_names_map[&row_genome];
        let col_idx = self.col_names_map[&col_genome];
        let idx = (row_idx as usize) * self.col_names.len() + (col_idx as usize);
        self.data[idx]
    }

    pub fn contains_matrix_by_names(&self, other: &Self) -> bool {
        let row_is_subset = other
            .row_names
            .iter()
            .all(|x| self.row_names_map.contains_key(x));
        let col_is_subset = other
            .col_names
            .iter()
            .all(|x| self.col_names_map.contains_key(x));
        row_is_subset && col_is_subset
    }

    pub fn update_from(&mut self, other: &Self) {
        assert!(self.contains_matrix_by_names(other));
        let (nrow, ncol) = other.shape();
        for i in 0..nrow {
            for j in 0..ncol {
                let v = other.get_by_positions(i as u32, j as u32);

                let row_idx = self.row_names_map[&other.row_names[i]];
                let col_idx = self.col_names_map[&other.col_names[j]];

                self.set_by_positions(row_idx, col_idx, v);
            }
        }
    }

    pub fn is_equal(&self, other: &Self) -> bool {
        (self.row_names == other.row_names)
            && (self.col_names == other.col_names)
            && (self.data == other.data)
    }

    pub fn transpose(&mut self) {
        let mut v = Vec::<T>::with_capacity(self.data.len());
        for col in 0..self.col_names.len() as u32 {
            for row in 0..self.row_names.len() as u32 {
                v.push(self.get_by_positions(row, col));
            }
        }
        self.data = v;
        std::mem::swap(&mut self.row_names, &mut self.col_names);
        std::mem::swap(&mut self.row_names_map, &mut self.col_names_map);
    }

    pub fn get_row_slice(&self, row_idx: u32) -> &[T] {
        let s = row_idx as usize * self.col_names.len();
        let e = (row_idx + 1) as usize * self.col_names.len();
        &self.data[s..e]
    }

    pub fn into_parts(mut self) -> (Vec<u32>, Vec<u32>, Vec<T>) {
        use std::mem::take;
        let row_names = take(&mut self.row_names);
        let col_names = take(&mut self.col_names);
        let data = take(&mut self.data);
        (row_names, col_names, data)
    }
}

impl<T> IntoParquet for NamedMatrix<T>
where
    T: Copy + Default + Sized + PartialEq,
    Vec<T>: IntoArrowArray,
    [T]: FromArrowArray,
{
    fn into_parquet(&mut self, p: impl AsRef<Path>) {
        use std::mem::take;
        // build array from genotype matrix
        let mat = take(&mut self.data).into_arrow_array();
        let row = take(&mut self.row_names).into_arrow_array();
        let col = take(&mut self.col_names).into_arrow_array();

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

impl<T> FromParquet for NamedMatrix<T>
where
    T: Copy + Default + Sized + PartialEq,
    Vec<T>: IntoArrowArray,
    [T]: FromArrowArray,
{
    fn from_parquet(p: impl AsRef<Path>) -> Self {
        let mut row_genomes = Vec::<u32>::new();
        let mut col_genomes = Vec::<u32>::new();
        let mut data = Vec::<T>::new();

        let file = File::open(&p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            let arr = record_batch.column(0).as_ref();
            let slice: &[T] = FromArrowArray::from_array_array(arr);
            data.extend_from_slice(slice);
        }

        let file = File::open(p.as_ref().with_extension("row")).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            row_genomes.extend_from_slice(FromArrowArray::from_array_array(record_batch.column(0)));
        }

        let file = File::open(p.as_ref().with_extension("col")).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut reader = builder.build().unwrap();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            col_genomes.extend_from_slice(FromArrowArray::from_array_array(record_batch.column(0)));
        }

        let mut matrix = Self::new(row_genomes, col_genomes);
        matrix.data = data;

        matrix
    }
}
