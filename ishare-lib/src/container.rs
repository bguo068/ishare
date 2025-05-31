pub mod histogram;
pub mod intervals;
pub mod intervaltree;

use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Histogram { source: histogram::Error },
    
    #[snafu(display("Intervals not sorted: intervals must be sorted before complement operation"))]
    IntervalsNotSorted,
    
    #[snafu(display("Invalid contig range: contig_min ({contig_min:?}) must be <= contig_max ({contig_max:?})"))]
    InvalidContigRange {
        contig_min: String,
        contig_max: String,
    },
    
    #[snafu(display("Interval out of bounds: interval [{start:?}, {end:?}) is outside contig range [{contig_min:?}, {contig_max:?})"))]
    IntervalOutOfBounds {
        start: String,
        end: String,
        contig_min: String,
        contig_max: String,
    },
}
