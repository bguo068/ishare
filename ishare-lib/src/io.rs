use arrow_array::{Array, ArrayRef, Float32Array, Float64Array, UInt32Array, UInt8Array};
use arrow_schema::ArrowError;
use snafu::prelude::*;
use std::backtrace::Backtrace;
use std::path::Path;
use std::sync::Arc;

type Result<T> = std::result::Result<T, Error>;

#[derive(Snafu, Debug)]
#[snafu(visibility(pub(crate)))]
pub enum Error {
    Downcast {
        backtrace: Option<Backtrace>,
    },
    StdIo {
        source: std::io::Error,
        backtrace: Option<Backtrace>,
    },
    Parquet {
        source: parquet::errors::ParquetError,
        backtrace: Option<Backtrace>,
    },
    Arrow {
        source: ArrowError,
        backtrace: Option<Backtrace>,
    },
}

pub trait IntoParquet {
    fn into_parquet(self, p: impl AsRef<Path>) -> Result<()>;
}

pub trait FromParquet: Sized {
    fn from_parquet(p: impl AsRef<Path>) -> Result<Self>;
}

pub trait IntoArrowArray {
    fn into_arrow_array(self) -> ArrayRef;
}

impl IntoArrowArray for Vec<u8> {
    fn into_arrow_array(self) -> ArrayRef {
        Arc::new(UInt8Array::from(self)) as ArrayRef
    }
}
impl IntoArrowArray for Vec<u32> {
    fn into_arrow_array(self) -> ArrayRef {
        Arc::new(UInt32Array::from(self)) as ArrayRef
    }
}
impl IntoArrowArray for Vec<f32> {
    fn into_arrow_array(self) -> ArrayRef {
        Arc::new(Float32Array::from(self)) as ArrayRef
    }
}
impl IntoArrowArray for Vec<f64> {
    fn into_arrow_array(self) -> ArrayRef {
        Arc::new(Float64Array::from(self)) as ArrayRef
    }
}

pub trait FromArrowArray {
    fn from_array_array(arr: &dyn Array) -> Result<&Self>;
}

impl FromArrowArray for [u8] {
    fn from_array_array(arr: &dyn Array) -> Result<&Self> {
        Ok(arr
            .as_any()
            .downcast_ref::<UInt8Array>()
            .context(DowncastSnafu {})?
            .values())
    }
}
impl FromArrowArray for [u32] {
    fn from_array_array(arr: &dyn Array) -> Result<&Self> {
        Ok(arr
            .as_any()
            .downcast_ref::<UInt32Array>()
            .context(DowncastSnafu {})?
            .values())
    }
}
impl FromArrowArray for [f32] {
    fn from_array_array(arr: &dyn Array) -> Result<&Self> {
        let x = arr
            .as_any()
            .downcast_ref::<Float32Array>()
            .context(DowncastSnafu {})?
            .values();
        Ok(x)
    }
}
impl FromArrowArray for [f64] {
    fn from_array_array(arr: &dyn Array) -> Result<&Self> {
        let x = arr
            .as_any()
            .downcast_ref::<Float64Array>()
            .context(DowncastSnafu {})?
            .values();
        Ok(x)
    }
}
