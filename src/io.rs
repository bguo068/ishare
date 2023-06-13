use std::path::Path;

pub trait IntoParquet {
    fn into_parquet(&mut self, p: impl AsRef<Path>);
}

pub trait FromParquet {
    fn from_parquet(p: impl AsRef<Path>) -> Self;
}
