use arrow_array::builder::GenericByteBuilder;
use arrow_array::types::GenericBinaryType;
use arrow_array::ArrayRef;
use arrow_array::GenericByteArray;
use arrow_array::RecordBatch;
use arrow_array::UInt32Array;
use arrow_schema::ArrowError;

use bstr::BString;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use snafu::prelude::*;
use std::backtrace::Backtrace;

type Result<T> = std::result::Result<T, Error>;

#[derive(Snafu, Debug)]
pub enum Error {
    Arrow {
        #[snafu(source(from(ArrowError, Box::new)))]
        source: Box<ArrowError>,
        backtrace: Box<Option<Backtrace>>,
    },
    StdIo {
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    Parquet {
        #[snafu(source(from(parquet::errors::ParquetError, Box::new)))]
        source: Box<parquet::errors::ParquetError>,
        backtrace: Box<Option<Backtrace>>,
    },
    BinarySearchNotFound {
        backtrace: Box<Option<Backtrace>>,
    },
    Downcast {
        backtrace: Box<Option<Backtrace>>,
    },
    MissingGenotypeData {
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Position {pos} is not greater than last position {last_pos} - sites must be added in order"))]
    PositionNotInOrder {
        pos: u32,
        last_pos: u32,
        backtrace: Box<Option<Backtrace>>,
    },
    #[snafu(display("Position {pos} is less than or equal to last position {last_pos} - sites must be added in strictly increasing order"))]
    PositionNotStrictlyIncreasing {
        pos: u32,
        last_pos: u32,
        backtrace: Box<Option<Backtrace>>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct Sites {
    gw_pos: Vec<u32>,    // gw_pos
    buf: Vec<u8>,        // alele string
    offsets: Vec<usize>, // offset of each allele starting byte
}
impl Sites {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn merge(&mut self, other: Self) {
        self.gw_pos.extend(other.gw_pos);
        let buf_start = self.buf.len();
        self.buf.extend(other.buf);
        self.offsets
            .extend(other.offsets.into_iter().map(|x| x + buf_start));
    }

    pub fn sort_by_position_then_allele(&mut self) -> Result<Vec<u32>> {
        let mut indices: Vec<u32> = (0..self.len() as u32).collect();

        indices.sort_by_key(|x| {
            (
                self.gw_pos[*x as usize],
                self.get_alleles_by_idx(*x as usize),
            )
        });

        let mut sites2 = Sites::new();

        for idx in indices.iter() {
            let (pos, bytes) = self.get_site_by_idx(*idx as usize);
            sites2.add_site_with_bytes(pos, bytes)?;
        }
        *self = sites2;

        Ok(indices)
    }

    /// without increasing the number of items
    pub fn append_bytes_to_last_allele(&mut self, bytes: &[u8]) {
        self.buf.extend_from_slice(bytes);
    }

    pub fn add_site_with_bytes(&mut self, pos: u32, bytes: &[u8]) -> Result<()> {
        self.offsets.push(self.buf.len());
        self.buf.extend_from_slice(bytes);
        //  site added in order
        if let Some(last_pos) = self.gw_pos.last() {
            ensure!(
                *last_pos <= pos,
                PositionNotInOrderSnafu {
                    pos,
                    last_pos: *last_pos
                }
            );
        }
        self.gw_pos.push(pos);
        Ok(())
    }

    pub fn add_site(&mut self, pos: u32, ab: &AlleleBuffer) -> Result<()> {
        self.offsets.push(self.buf.len());
        let mut first_allele = true;
        for i in 0..ab.len() {
            let bytes = ab.get(i);
            let enc = ab.get_enc(i as u8);
            // skip encs that are marked for cleaning up (ALT with AC=0)
            if enc == u8::MAX {
                continue;
            }
            // Add space separator before non-first alleles
            if !first_allele {
                self.buf.extend_from_slice(b" ");
            }
            self.buf.extend_from_slice(bytes);
            first_allele = false;
        }
        // ensure site added in order
        if let Some(last_pos) = self.gw_pos.last() {
            ensure!(
                *last_pos < pos,
                PositionNotStrictlyIncreasingSnafu {
                    pos,
                    last_pos: *last_pos
                }
            );
        }
        self.gw_pos.push(pos);
        Ok(())
    }

    /// When positions are known to be sorted and unique
    pub fn get_site_by_position(&self, pos: u32) -> Result<&[u8]> {
        let idx = self
            .gw_pos
            .binary_search(&pos)
            .map_err(|_| BinarySearchNotFoundSnafu {}.build())?;
        let s = self.offsets[idx];
        let mut e = self.buf.len();
        if idx < self.offsets.len() - 1 {
            e = self.offsets[idx + 1];
        }
        Ok(&self.buf[s..e])
    }

    pub fn get_gw_pos_slice(&self) -> &[u32] {
        &self.gw_pos[..]
    }

    // When positions are sorted but may not be unique
    pub fn get_idx_by_position(&self, pos: u32) -> (usize, usize) {
        (
            self.gw_pos.partition_point(|x| *x < pos),
            self.gw_pos.partition_point(|x| *x <= pos),
        )
    }

    pub fn get_alleles_by_idx(&self, idx: usize) -> &[u8] {
        let s = self.offsets[idx];
        let mut e = self.buf.len();
        if idx < self.offsets.len() - 1 {
            e = self.offsets[idx + 1];
        }
        &self.buf[s..e]
    }

    pub fn get_site_by_idx(&self, idx: usize) -> (u32, &[u8]) {
        let s = self.offsets[idx];
        let mut e = self.buf.len();
        if idx < self.offsets.len() - 1 {
            e = self.offsets[idx + 1];
        }
        // println!("get_site_by_idx: s={s},e={e}, buflen={}", self.buf.len());
        (self.gw_pos[idx], &self.buf[s..e])
    }

    pub fn len(&self) -> usize {
        self.gw_pos.len()
    }

    pub fn is_empty(&self) -> bool {
        self.gw_pos.len() == 0
    }

    pub fn into_parquet_file(self, p: impl AsRef<Path>) -> Result<()> {
        // for col2
        let mut builder = GenericByteBuilder::<GenericBinaryType<i32>>::new();
        for i in 0..self.len() {
            let (_, s) = self.get_site_by_idx(i);
            builder.append_value(s);
        }
        // array
        let col1 = UInt32Array::from(self.gw_pos);
        let col2 = builder.finish();
        // record batch
        let batch = RecordBatch::try_from_iter(vec![
            ("gw_pos", Arc::new(col1) as ArrayRef),
            ("alleles", Arc::new(col2) as ArrayRef),
        ])
        .context(ArrowSnafu {})?;
        // writer
        let file = File::create(p.as_ref()).context(StdIoSnafu {})?;
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
        let file = File::open(p).context(StdIoSnafu {})?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).context(ParquetSnafu {})?;
        let mut reader = builder.build().context(ParquetSnafu {})?;

        let mut sites = Sites::new();

        for record_batch in &mut reader {
            let record_batch = record_batch.context(ArrowSnafu {})?;
            let pos_iter = record_batch
                .column(0)
                .as_any()
                .downcast_ref::<UInt32Array>()
                .context(DowncastSnafu {})?
                .values()
                .iter()
                .copied();
            let bytes_iter = record_batch
                .column(1)
                .as_any()
                .downcast_ref::<GenericByteArray<GenericBinaryType<i32>>>()
                .context(DowncastSnafu {})?
                .into_iter();

            for (pos, bytes) in pos_iter.zip(bytes_iter) {
                sites.add_site_with_bytes(pos, bytes.context(MissingGenotypeDataSnafu {})?)?;
            }
        }

        Ok(sites)
    }
}

#[derive(Debug, Default)]
pub struct AlleleBuffer {
    data: BString,
    offsets: Vec<usize>,
    enc: Vec<u8>,
}
impl AlleleBuffer {
    pub fn new() -> Self {
        Self {
            data: BString::new(vec![0; 32]),
            offsets: Vec::new(),
            enc: Vec::new(),
        }
    }
    pub fn push(&mut self, allele: &[u8]) {
        self.enc.push(self.offsets.len() as u8);
        self.offsets.push(self.data.len());
        self.data.extend_from_slice(allele);
    }
    pub fn push_to_data_only(&mut self, allele: &[u8]) {
        self.data.extend_from_slice(allele);
    }
    pub fn get(&self, i: usize) -> &[u8] {
        let s = self.offsets[i];
        let mut e = self.data.len();
        if i < self.offsets.len() - 1 {
            e = self.offsets[i + 1];
        }
        &self.data[s..e]
    }
    pub fn get_enc(&self, allele_ix_in_vcf: u8) -> u8 {
        self.enc[allele_ix_in_vcf as usize]
    }
    pub fn set_enc(&mut self, allele_ix_in_vcf: u8, new_ix: u8) {
        self.enc[allele_ix_in_vcf as usize] = new_ix;
    }
    pub fn len(&self) -> usize {
        self.offsets.len()
    }
    pub fn is_empty(&self) -> bool {
        self.offsets.len() == 0
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.offsets.clear();
        self.enc.clear();
    }
}
