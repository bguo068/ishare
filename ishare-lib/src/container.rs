pub mod histogram;
pub mod intervals;
pub mod intervaltree;

use std::backtrace::Backtrace;

use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Histogram {
        // non leaf
        #[snafu(backtrace)]
        source: histogram::Error,
    },

    #[snafu(display(
        "Intervals not sorted: intervals must be sorted before complement operation"
    ))]
    IntervalsNotSorted {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },

    #[snafu(display("Invalid contig range: {msg})"))]
    InvalidContigRange {
        // leaf
        msg: Box<String>,
        backtrace: Box<Option<Backtrace>>,
    },

    #[snafu(display("Interval out of bounds: {msg}"))]
    IntervalOutOfBounds {
        // leaf
        msg: Box<String>,
        backtrace: Box<Option<Backtrace>>,
    },
}
