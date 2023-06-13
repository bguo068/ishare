use std::cmp::{PartialEq, PartialOrd};

pub struct Histogram<T>
where
    T: PartialEq + PartialOrd + Copy,
{
    bins: Vec<T>,
    counts: Vec<usize>,
}

impl<'a, T> Histogram<T>
where
    T: PartialEq + PartialOrd + Copy + 'a,
{
    pub fn new(bin_iter: impl Iterator<Item = &'a T>) -> Self {
        let mut bins: Vec<T> = bin_iter.map(|x| *x).collect();
        bins.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            bins.iter().zip(bins.iter().skip(1)).all(|(a, b)| *a < *b),
            "boundaries should be sorted and unique"
        );
        assert!(bins.len() >= 1, "min length of boundaries is 1");
        let n = bins.len();
        Self {
            bins,
            counts: vec![0; n],
        }
    }

    pub fn analyze(&mut self, val_iter: impl Iterator<Item = &'a T>) {
        let min = self.bins[0];
        for val in val_iter {
            if val.partial_cmp(&min).unwrap().is_lt() {
                // ignore those that are two short/small
                continue;
            }
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
    let val = vec![0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4];
    let mut hist = Histogram::new(bins.iter());
    hist.analyze(val.iter());

    assert_eq!(hist.get_counts(), &vec![2, 3, 9]);
}
