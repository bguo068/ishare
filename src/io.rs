use arrow::array::*;
use std::path::Path;
use std::sync::Arc;

pub trait IntoParquet {
    fn into_parquet(&mut self, p: impl AsRef<Path>);
}

pub trait FromParquet {
    fn from_parquet(p: impl AsRef<Path>) -> Self;
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
    fn from_array_array<'a>(arr: &'a dyn Array) -> &'a Self;
}

impl FromArrowArray for [u8] {
    fn from_array_array<'a>(arr: &'a dyn Array) -> &'a Self {
        arr.as_any().downcast_ref::<UInt8Array>().unwrap().values()
    }
}
impl FromArrowArray for [u32] {
    fn from_array_array<'a>(arr: &'a dyn Array) -> &'a Self {
        arr.as_any().downcast_ref::<UInt32Array>().unwrap().values()
    }
}
impl FromArrowArray for [f32] {
    fn from_array_array<'a>(arr: &'a dyn Array) -> &'a Self {
        let x = arr
            .as_any()
            .downcast_ref::<Float32Array>()
            .unwrap()
            .values();
        x
    }
}
impl FromArrowArray for [f64] {
    fn from_array_array<'a>(arr: &'a dyn Array) -> &'a Self {
        let x = arr
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap()
            .values();
        x
    }
}
