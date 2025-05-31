pub mod xirs;

use snafu::prelude::*;

#[derive(Snafu, Debug)]
pub enum Error {
    #[snafu(transparent)]
    Xirs { source: xirs::Error },
}
