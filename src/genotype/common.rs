use arrow::array::*;
use arrow::buffer::BooleanBuffer;

use crate::site::Sites;
use arrow::record_batch::RecordBatch;
use bitvec::prelude::*;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

// Matrix of 0 and 1.
///
/// Multiallelic sites must be presented as multiple lines of biallelic genotypes
#[derive(Serialize, Deserialize)]
pub struct GenotypeMatrix {
    bv: BitVec<u64, Lsb0>,
    ncols: usize,
}

impl GenotypeMatrix {
    pub fn new(ncol: usize) -> Self {
        Self {
            bv: BitVec::<u64, Lsb0>::new(),
            ncols: ncol,
        }
    }

    pub fn nrows(&self) -> usize {
        let nrow = self.bv.len() / self.ncols;
        assert_eq!(self.ncols * nrow, self.bv.len());
        nrow
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
        for col in 0..self.nrows() {
            for row in 0..self.ncols() {
                bv2.push(self.get_at(row, col))
            }
        }
        let ncols2 = self.nrows();
        Self {
            bv: bv2,
            ncols: ncols2,
        }
    }

    pub fn append_row(&mut self, other: &BitSlice<u64>) {
        let len1 = self.bv.len();
        self.bv.extend_from_bitslice(other);
        let len2 = self.bv.len();
        assert!(len1 + self.ncols == len2);
    }

    pub fn reorder_rows(self, row_orders: &[u32]) -> Self {
        assert_eq!(row_orders.len(), self.nrows());
        let mut gt2 = GenotypeMatrix::new(self.ncols);
        gt2.bv.reserve(self.bv.len());

        for idx in row_orders {
            let row_slice = self.get_row(*idx as usize);
            gt2.append_row(row_slice);
        }
        gt2
    }

    pub fn merge(&mut self, other: Self) {
        assert_eq!(self.ncols, other.ncols);
        let slice = other.bv.as_bitslice();
        self.bv.extend_from_bitslice(slice);
    }

    pub fn has_allele(&self, pos: u32, gid: u32, allele: &[u8], sites: &Sites) -> bool {
        let (s, e) = sites.get_idx_by_position(pos);
        if s == e {
            false
        } else {
            match allele {
                b"REF" => (s..e).all(|pos_idx| !self.get_at(pos_idx, gid as usize)),
                _ => (s..e).any(|pos_idx| {
                    (sites.get_alleles_by_idx(pos_idx) == allele)
                        && (self.get_at(pos_idx, gid as usize))
                }),
            }
        }
    }
    pub fn into_parquet_file(self, p: impl AsRef<Path>) {
        // build array from genotype matrix
        let mut builder = BooleanBufferBuilder::new(self.bv.len());
        for i in 0..self.bv.len() {
            builder.append(self.bv[i]);
        }
        let buffer: BooleanBuffer = builder.finish();
        let arr = BooleanArray::new(buffer, None);

        // trick: use field name to save the ncol member variable
        let fieldname = format!("{}", self.ncols);

        let batch =
            RecordBatch::try_from_iter(vec![(fieldname, Arc::new(arr) as ArrayRef)]).unwrap();
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
        let ncol = builder.schema().fields[0].name().parse::<usize>().unwrap();
        let mut reader = builder.build().unwrap();
        let mut gm = GenotypeMatrix::new(ncol);
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            record_batch
                .column(0)
                .as_any()
                .downcast_ref::<BooleanArray>()
                .unwrap()
                .into_iter()
                .for_each(|x| gm.bv.push(x.unwrap()));
        }

        gm
    }

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
}
