pub mod histogram;
pub mod intervals;
pub mod intervaltree;

use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Histogram { source: histogram::Error },

    #[snafu(display(
        "Intervals not sorted: intervals must be sorted before complement operation"
    ))]
    IntervalsNotSorted,
    #[snafu(display("Invalid contig range: {msg})"))]
    InvalidContigRange { msg: Box<String> },

    #[snafu(display("Interval out of bounds: {msg}"))]
    IntervalOutOfBounds { msg: Box<String> },
}
