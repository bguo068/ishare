use super::super::traits::TotalOrd;
use snafu::{ensure, Snafu};
use std::backtrace::Backtrace;

#[derive(Debug, Snafu)]
pub enum Error {
    // InComparableError { backtrace: Option<Backtrace> },
    BoundaryNotUniqueError { backtrace: Option<Backtrace> },
    EmptyBoundaryError { backtrace: Option<Backtrace> },
}
type Result<T> = std::result::Result<T, Error>;

pub struct Histogram<T>
where
    T: TotalOrd + Copy,
{
    bins: Vec<T>,
    counts: Vec<usize>,
}

impl<'a, T> Histogram<T>
where
    T: TotalOrd + Copy + 'a,
{
    pub fn new(bin_iter: impl Iterator<Item = &'a T>) -> Result<Self> {
        let mut bins: Vec<T> = bin_iter.copied().collect();

        // boundary should not be empty
        ensure!(!bins.is_empty(), EmptyBoundarySnafu {});

        bins.sort_by(|a, b| a.total_cmp(b));

        ensure!(
            bins.iter().zip(bins.iter().skip(1)).all(|(a, b)| *a != *b),
            BoundaryNotUniqueSnafu {}
        );

        let n = bins.len();
        Ok(Self {
            bins,
            counts: vec![0; n],
        })
    }

    pub fn analyze(&mut self, val_iter: impl Iterator<Item = &'a T>) {
        let min = self.bins[0];
        for val in val_iter {
            if val.total_cmp(&min).is_lt() {
                // ignore those that are two short/small
                continue;
            }
            // unwrap won't panic as both val_iter and boundary both have
            // been tested for comparability
            let idx = self
                .bins
                .partition_point(|x| x.partial_cmp(val).unwrap().is_le());
            self.counts[idx - 1] += 1;
        }
    }

    pub fn get_bins(&self) -> &Vec<T> {
        &self.bins
    }
    pub fn get_counts(&self) -> &Vec<usize> {
        &self.counts
    }
}

#[test]
fn test_histogram() {
    let bins = [1, 2, 3];
    let val = [0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4];
    let mut hist = Histogram::new(bins.iter()).unwrap();
    hist.analyze(val.iter());

    assert_eq!(hist.get_counts(), &vec![2, 3, 9]);
}
