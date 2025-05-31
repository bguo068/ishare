#![cfg_attr(not(test), warn(clippy::unwrap_used))]
#![cfg_attr(not(test), warn(clippy::expect_used))]

pub mod container;
pub mod genome;
pub mod genotype;
pub mod gmap;
pub mod indiv;
pub mod io;
pub mod rfmix;
pub mod share;
pub mod site;
pub mod stat;
#[cfg(test)]
pub mod tests;
pub mod traits;
pub mod utils;
pub mod vcf;
