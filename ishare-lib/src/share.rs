pub mod ibd;
pub mod mat;
pub mod pairs;

use snafu::prelude::*;

#[derive(Snafu, Debug)]
pub enum Error {
    #[snafu(transparent)]
    Mat { source: mat::Error },
    #[snafu(transparent)]
    Ibd { source: ibd::Error },
}
