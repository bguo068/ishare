pub mod afreq;
pub mod common;
pub mod rare;

use snafu::prelude::*;

#[derive(Debug, Snafu)]
pub enum Error {
    #[snafu(transparent)]
    Afreq { source: afreq::Error },
    #[snafu(transparent)]
    Common { source: common::Error },
    #[snafu(transparent)]
    Rare { source: rare::Error },
}
