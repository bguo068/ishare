pub mod coverage;
pub mod ibdseg;
pub mod ibdset;
pub mod overlap;
pub mod peak;

use snafu::prelude::*;

#[derive(Snafu, Debug)]
pub enum Error {
    #[snafu(display("Failed to read directory: {}", path.display()))]
    ReadDirectory {
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
    },

    #[snafu(display("Failed to create file: {}", path.display()))]
    CreateFile {
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
    },

    #[snafu(display("Failed to write to file: {}", path.display()))]
    WriteFile {
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
    },

    #[snafu(display("Invalid filename: {}", filename))]
    InvalidFilename { filename: Box<String> },

    #[snafu(display("Failed to parse string as UTF-8"))]
    Utf8Parse { source: std::str::Utf8Error },

    #[snafu(display("Failed to parse: {msg}"))]
    //'{value}' as {type_name}"
    ParseValue {
        source: std::num::ParseFloatError,
        msg: Box<String>,
    },

    #[snafu(display("Failed to parse '{value}' as integer"))]
    ParseInt {
        source: std::num::ParseIntError,
        value: Box<String>,
    },

    #[snafu(display("Sample '{sample}' not found in individuals"))]
    SampleNotFound { sample: Box<String> },

    #[snafu(display("Invalid haplotype value: {value}"))]
    InvalidHaplotype { value: u8 },

    #[snafu(display("CSV reading error"))]
    CsvRead { source: csv::Error },

    #[snafu(display("BGZ reading error"))]
    BgzRead {
        source: Box<rust_htslib::errors::Error>,
    },

    #[snafu(display("Thread pool creation failed"))]
    ThreadPool {
        source: Box<rust_htslib::errors::Error>,
    },

    #[snafu(display("Array index out of bounds"))]
    IndexOutOfBounds,

    #[snafu(display("Missing required data"))]
    MissingData,

    #[snafu(display("Invalid enum value: {value}"))]
    InvalidEnumValue { value: Box<String> },

    #[snafu(display("Container operation failed: {details}"))]
    ContainerOperation { details: Box<String> },
}
