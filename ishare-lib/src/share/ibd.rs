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
        path: std::path::PathBuf,
    },

    #[snafu(display("Failed to create file: {}", path.display()))]
    CreateFile {
        source: std::io::Error,
        path: std::path::PathBuf,
    },

    #[snafu(display("Failed to write to file: {}", path.display()))]
    WriteFile {
        source: std::io::Error,
        path: std::path::PathBuf,
    },

    #[snafu(display("Invalid filename: {}", filename))]
    InvalidFilename { filename: String },

    #[snafu(display("Failed to parse string as UTF-8"))]
    Utf8Parse { source: std::str::Utf8Error },

    #[snafu(display("Failed to parse '{value}' as {type_name}"))]
    ParseValue {
        source: std::num::ParseFloatError,
        value: String,
        type_name: String,
    },

    #[snafu(display("Failed to parse '{value}' as integer"))]
    ParseInt {
        source: std::num::ParseIntError,
        value: String,
    },

    #[snafu(display("Sample '{sample}' not found in individuals"))]
    SampleNotFound { sample: String },

    #[snafu(display("Invalid haplotype value: {value}"))]
    InvalidHaplotype { value: u8 },

    #[snafu(display("CSV reading error"))]
    CsvRead { source: csv::Error },

    #[snafu(display("BGZ reading error"))]
    BgzRead { source: rust_htslib::errors::Error },

    #[snafu(display("Thread pool creation failed"))]
    ThreadPool { source: rust_htslib::errors::Error },

    #[snafu(display("Array index out of bounds"))]
    IndexOutOfBounds,

    #[snafu(display("Missing required data"))]
    MissingData,

    #[snafu(display("Invalid enum value: {value}"))]
    InvalidEnumValue { value: String },

    #[snafu(display("Container operation failed: {details}"))]
    ContainerOperation { details: String },
}
