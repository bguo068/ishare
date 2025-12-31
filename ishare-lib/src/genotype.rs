pub mod afreq;
pub mod common;
pub mod rare;

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    Afreq {
        // non leaf
        #[snafu(backtrace)]
        source: afreq::Error,
    },
    // #[snafu(transparent)]
    Common {
        // non leaf
        #[snafu(backtrace)]
        source: common::Error,
    },
    // #[snafu(transparent)]
    Rare {
        // non leaf
        #[snafu(backtrace)]
        source: rare::Error,
    },
}
