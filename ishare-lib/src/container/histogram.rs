use super::super::traits::TotalOrd;
use snafu::{ensure, Snafu};
use std::backtrace::Backtrace;

#[derive(Debug, Snafu)]
pub enum Error {
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
            let idx = self.bins.partition_point(|x| x.total_cmp(val).is_le());
            self.counts[idx - 1] += 1;
        }
    }

    pub fn get_bins(&self) -> &[T] {
        &self.bins
    }
    pub fn get_counts(&self) -> &[usize] {
        &self.counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_histogram_basic() {
        let bins = [1, 2, 3];
        let val = [0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4];
        let mut hist = Histogram::new(bins.iter()).unwrap();
        hist.analyze(val.iter());

        assert_eq!(hist.get_counts(), &vec![2, 3, 9]);
    }

    mod boundary_validation {
        use super::*;

        #[test]
        fn test_empty_boundary_error() {
            let empty_bins: Vec<i32> = vec![];
            let result = Histogram::new(empty_bins.iter());
            assert!(matches!(result, Err(Error::EmptyBoundaryError { .. })));
        }

        #[test]
        fn test_single_boundary() {
            let bins = [42];
            let val = [41, 42, 43];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Single boundary creates one bin, values >= boundary go in that bin
            assert_eq!(hist.get_counts(), &vec![2]); // 42, 43
        }

        #[test]
        fn test_duplicate_boundaries_error() {
            let bins = [1, 2, 2, 3];
            let result = Histogram::new(bins.iter());
            assert!(matches!(result, Err(Error::BoundaryNotUniqueError { .. })));
        }

        #[test]
        fn test_unordered_boundaries() {
            let bins = [3, 1, 2];
            let val = [0, 1, 2, 3, 4];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Should be sorted internally: [1, 2, 3]
            assert_eq!(hist.get_bins(), &vec![1, 2, 3]);
            assert_eq!(hist.get_counts(), &vec![1, 1, 2]); // [1], [2], [3,4]
        }

        #[test]
        fn test_all_identical_boundaries_error() {
            let bins = [5, 5, 5];
            let result = Histogram::new(bins.iter());
            assert!(matches!(result, Err(Error::BoundaryNotUniqueError { .. })));
        }
    }

    mod edge_case_analysis {
        use super::*;

        #[test]
        fn test_values_below_minimum() {
            let bins = [10, 20, 30];
            let val = [1, 2, 5, 15, 25, 35];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Values below minimum (1,2,5) are ignored
            assert_eq!(hist.get_counts(), &vec![1, 1, 1]); // [15], [25], [35]
        }

        #[test]
        fn test_values_on_boundaries() {
            let bins = [10, 20, 30];
            let val = [10, 20, 30];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Values exactly on boundaries
            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }

        #[test]
        fn test_values_above_maximum() {
            let bins = [10, 20, 30];
            let val = [15, 25, 35, 40, 50];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // All values above maximum go in last bin
            assert_eq!(hist.get_counts(), &vec![1, 1, 3]); // [15], [25], [35,40,50]
        }

        #[test]
        fn test_empty_values() {
            let bins = [1, 2, 3];
            let val: Vec<i32> = vec![];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![0, 0, 0]);
        }

        #[test]
        fn test_all_values_below_minimum() {
            let bins = [100, 200, 300];
            let val = [1, 2, 3, 50, 99];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![0, 0, 0]);
        }

        #[test]
        fn test_large_dataset() {
            let bins = [0, 25, 50, 75, 100];
            let val: Vec<i32> = (0..1000).collect();
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Values: 0-999
            // bin 0: [0-24] = 25 values (0 to 24)
            // bin 1: [25-49] = 25 values (25 to 49)
            // bin 2: [50-74] = 25 values (50 to 74)
            // bin 3: [75-99] = 25 values (75 to 99)
            // bin 4: [100-999] = 900 values (100 to 999)
            assert_eq!(hist.get_counts(), &vec![25, 25, 25, 25, 900]);
        }
    }

    mod different_data_types {
        use super::*;

        #[test]
        fn test_u32_histogram() {
            let bins: Vec<u32> = vec![10, 20, 30];
            let val: Vec<u32> = vec![5, 15, 25, 35];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]); // [15], [25], [35]
        }

        #[test]
        fn test_u64_histogram() {
            let bins: Vec<u64> = vec![1000000, 2000000, 3000000];
            let val: Vec<u64> = vec![500000, 1500000, 2500000, 3500000];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }

        #[test]
        fn test_f32_histogram() {
            let bins: Vec<f32> = vec![1.0, 2.0, 3.0];
            let val: Vec<f32> = vec![0.5, 1.5, 2.5, 3.5];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }

        #[test]
        fn test_f64_histogram() {
            let bins: Vec<f64> = vec![0.1, 0.2, 0.3];
            let val: Vec<f64> = vec![0.05, 0.15, 0.25, 0.35];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }

        #[test]
        fn test_mixed_precision_values() {
            let bins: Vec<f64> = vec![1.1, 2.2, 3.3];
            let val: Vec<f64> = vec![1.05, 1.15, 2.15, 2.25, 3.25, 3.35];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // 1.05: below minimum, ignored
            // 1.15: bin 0 (>= 1.1, < 2.2)
            // 2.15: bin 0 (>= 1.1, < 2.2)
            // 2.25: bin 1 (>= 2.2, < 3.3)
            // 3.25: bin 1 (>= 2.2, < 3.3)
            // 3.35: bin 2 (>= 3.3)
            assert_eq!(hist.get_counts(), &vec![2, 2, 1]);
        }

        #[test]
        fn test_tuple_histogram() {
            let bins: Vec<(u32, u32)> = vec![(1, 0), (2, 0), (3, 0)];
            let val: Vec<(u32, u32)> = vec![(0, 5), (1, 5), (2, 5), (3, 5), (4, 5)];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // Tuples ordered lexicographically
            assert_eq!(hist.get_counts(), &vec![1, 1, 2]);
        }
    }

    mod boundary_conditions {
        use super::*;

        #[test]
        fn test_zero_length_intervals() {
            let bins = [0, 0]; // This should fail due to duplicate boundary
            let result = Histogram::new(bins.iter());
            assert!(matches!(result, Err(Error::BoundaryNotUniqueError { .. })));
        }

        #[test]
        fn test_maximum_value_boundaries() {
            let bins = [u32::MAX - 2, u32::MAX - 1, u32::MAX];
            let val = [u32::MAX - 3, u32::MAX - 1, u32::MAX];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // u32::MAX - 3: below minimum, ignored
            // u32::MAX - 1: falls in bin 1 (>= u32::MAX - 1, < u32::MAX)
            // u32::MAX: falls in bin 2 (>= u32::MAX)
            assert_eq!(hist.get_counts(), &vec![0, 1, 1]);
        }

        #[test]
        fn test_minimum_value_boundaries() {
            let bins = [0u32, 1, 2];
            let val = [0, 1, 2];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }

        #[test]
        fn test_negative_values() {
            let bins = [-10, 0, 10];
            let val = [-20, -5, 5, 15];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]); // [-5], [5], [15]
        }

        #[test]
        fn test_very_close_boundaries() {
            let bins = [1.0, 1.0000001, 1.0000002];
            let val = [1.00000005, 1.00000015, 1.00000025];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            assert_eq!(hist.get_counts(), &vec![1, 1, 1]);
        }
    }

    mod comprehensive_scenarios {
        use super::*;

        #[test]
        fn test_statistical_distribution() {
            // Simulate normal distribution-like data
            let bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
            let val = [
                5, 15, 25, 35, 45, 55, 65, 75, 85, 95, // One per bin
                25, 35, 45, 55, 65, // Extra in middle bins
                45, 55, // Even more in center
            ];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // With 11 bins, we have 10 ranges: [0-9], [10-19], [20-29], [30-39], [40-49], [50-59], [60-69], [70-79], [80-89], [90-99], [100+]
            // 5: bin 0, 15: bin 1, 25,25: bin 2, 35,35: bin 3, 45,45,45: bin 4, 55,55,55: bin 5, 65,65: bin 6, 75: bin 7, 85: bin 8, 95: bin 9
            let counts = hist.get_counts();
            assert_eq!(counts.len(), 11); // 11 bins from 11 boundaries
            assert_eq!(counts[0], 1); // [5]
            assert_eq!(counts[1], 1); // [15]
            assert_eq!(counts[2], 2); // [25, 25]
            assert_eq!(counts[3], 2); // [35, 35]
            assert_eq!(counts[4], 3); // [45, 45, 45]
            assert_eq!(counts[5], 3); // [55, 55, 55]
            assert_eq!(counts[6], 2); // [65, 65]
            assert_eq!(counts[7], 1); // [75]
            assert_eq!(counts[8], 1); // [85]
            assert_eq!(counts[9], 1); // [95]
            assert_eq!(counts[10], 0); // [>=100] - none in this range
        }

        #[test]
        fn test_boundary_edge_cases_comprehensive() {
            let bins = [1, 5, 10];
            let val = [0, 1, 2, 4, 5, 6, 9, 10, 11, 15];
            let mut hist = Histogram::new(bins.iter()).unwrap();
            hist.analyze(val.iter());

            // 0: ignored (below min)
            // 1,2,4: bin 0 (>= 1, < 5)
            // 5,6,9: bin 1 (>= 5, < 10)
            // 10,11,15: bin 2 (>= 10)
            assert_eq!(hist.get_counts(), &vec![3, 3, 3]);
        }
    }
}
