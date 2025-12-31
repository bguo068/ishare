use arrow_array::{ArrayRef, Float32Array, RecordBatch};
use arrow_schema::ArrowError;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::backtrace::Backtrace;
use std::fs::File;
use std::sync::Arc;

use parquet::errors::ParquetError;
use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    Parquet {
        // non leaf
        #[snafu(source(from(ParquetError, Box::new)))]
        source: Box<ParquetError>,
        backtrace: Box<Option<Backtrace>>,
    },
    // #[snafu(transparent)]
    Arrow {
        // non leaf
        #[snafu(source(from(ArrowError, Box::new)))]
        source: Box<ArrowError>,
        backtrace: Box<Option<Backtrace>>,
    },
    StdIo {
        // non leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
}
type Result<T> = std::result::Result<T, Error>;

pub fn write_pair_total(
    total1_vec: Vec<f32>,
    total2_vec: Vec<f32>,
    colname1: &str,
    colname2: &str,
    p: impl AsRef<std::path::Path>,
) -> Result<()> {
    let t1 = Float32Array::from(total1_vec);
    let t2 = Float32Array::from(total2_vec);
    let batch = RecordBatch::try_from_iter(vec![
        (colname1.to_owned(), Arc::new(t1) as ArrayRef),
        (colname2.to_owned(), Arc::new(t2) as ArrayRef),
    ])
    .context(ArrowSnafu)?;

    let file = File::create(p.as_ref()).context(StdIoSnafu)?;

    // Default writer properties
    let props = WriterProperties::builder().build();

    let mut writer =
        ArrowWriter::try_new(file, batch.schema(), Some(props)).context(ParquetSnafu)?;

    writer.write(&batch).context(ParquetSnafu)?;

    // writer must be closed to write footer
    writer.close().context(ParquetSnafu)?;
    Ok(())
}
