pub mod xirs;

use snafu::prelude::*;

#[derive(Snafu, Debug)]
pub enum Error {
    #[snafu(transparent)]
    Xirs {
        // non-leaf
        #[snafu(backtrace)]
        source: xirs::Error,
    },
}
