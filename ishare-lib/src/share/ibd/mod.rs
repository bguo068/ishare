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
        // leaf
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Failed to create file: {}", path.display()))]
    CreateFile {
        // leaf
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Failed to write to file: {}", path.display()))]
    WriteFile {
        // leaf
        source: std::io::Error,
        path: Box<std::path::PathBuf>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Invalid filename: {}", filename))]
    InvalidFilename {
        // leaf
        filename: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Failed to parse string as UTF-8"))]
    Utf8Parse {
        // leaf
        source: std::str::Utf8Error,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Failed to parse: {msg}"))]
    ParseValue {
        // leaf
        source: std::num::ParseFloatError,
        msg: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Failed to parse '{value}' as integer"))]
    ParseInt {
        // leaf
        source: std::num::ParseIntError,
        value: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Sample '{sample}' not found in individuals"))]
    SampleNotFound {
        // leaf
        sample: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Invalid haplotype value: {value}"))]
    InvalidHaplotype {
        // leaf
        value: u8,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("CSV reading error"))]
    CsvRead {
        // leaf
        source: csv::Error,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("BGZ reading error"))]
    BgzRead {
        // leaf
        source: Box<rust_htslib::errors::Error>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Thread pool creation failed"))]
    ThreadPool {
        // leaf
        source: Box<rust_htslib::errors::Error>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Array index out of bounds"))]
    IndexOutOfBounds {
        // leaf
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Missing required data"))]
    MissingData {
        // leaf
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Invalid enum value: {value}"))]
    InvalidEnumValue {
        // leaf
        value: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },

    #[snafu(display("Container operation failed: {details}"))]
    ContainerOperation {
        // leaf
        details: Box<String>,
        backtrace: Box<std::backtrace::Backtrace>,
    },
}
