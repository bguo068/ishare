pub mod ibd;
pub mod mat;
pub mod pairs;

use snafu::prelude::*;

#[derive(Snafu, Debug)]
pub enum Error {
    #[snafu(transparent)]
    Mat {
        // non leaf
        #[snafu(backtrace)]
        source: mat::Error,
    },
    #[snafu(transparent)]
    Ibd {
        // non leaf
        #[snafu(backtrace)]
        source: ibd::Error,
    },
}
