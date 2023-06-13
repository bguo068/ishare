use std::fmt::Debug;
use std::ops::Range;

#[derive(Clone, Debug)]
pub struct Intervals<T: PartialOrd<T> + Copy + Default + Debug>(Vec<Range<T>>);

impl<'a, T: PartialOrd<T> + Copy + Default + Debug> Intervals<T> {
    pub fn interval_is_less(a: &Range<T>, b: &Range<T>) -> bool {
        (a.start, a.end) < (b.start, b.end)
    }
    pub fn interval_is_equal(a: &Range<T>, b: &Range<T>) -> bool {
        (a.start, a.end) == (b.start, b.end)
    }
    pub fn is_equal(&self, other: &Self) -> bool {
        if self.0.len() != other.0.len() {
            false
        } else {
            self.0
                .iter()
                .zip(other.0.iter())
                .all(|(a, b)| Self::interval_is_equal(a, b))
        }
    }
    pub fn new() -> Self {
        Intervals(Vec::new())
    }
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn push(&mut self, r: Range<T>) {
        self.0.push(r);
    }

    pub fn extend_from_iter(&mut self, t: impl Iterator<Item = (T, T)>) {
        for e in t {
            self.push(e.0..e.1);
        }
    }

    pub fn from_tuples(t: &[(T, T)]) -> Self {
        let mut intervals = Intervals::new();
        for e in t {
            intervals.push(e.0..e.1);
        }
        intervals
    }
    pub fn sort(&mut self) {
        self.0
            .sort_by(|a, b| (a.start, a.end).partial_cmp(&(b.start, b.end)).unwrap());
    }
    pub fn is_sorted(&self) -> bool {
        self.0
            .iter()
            .zip(self.0.iter().skip(1))
            .all(|(a, b)| Self::interval_is_less(a, b))
    }
    pub fn merge(&mut self) {
        if self.0.len() < 2 {
            return;
        }
        if !self.is_sorted() {
            self.sort();
        }
        let (mut e, mut rest) = self.0.split_first_mut().unwrap();
        while rest.len() > 0 {
            if e.end >= rest[0].start {
                rest[0].start = e.start;
                e.end = e.start;
            }
            (e, rest) = rest.split_first_mut().unwrap();
        }
        self.0.retain(|r| r.start != r.end);
    }

    pub fn complement(&mut self, min: T, max: T) {
        self.merge();
        if self.0.len() == 0 {
            self.push(min..max);
        } else {
            let the_min = self.0.first().unwrap().start;
            let the_max = self.0.last().unwrap().end;
            assert!(the_min >= min);
            assert!(the_max <= max);

            let (mut e, mut rest) = self.0.split_first_mut().unwrap();
            while rest.len() > 0 {
                e.start = e.end;
                e.end = rest[0].start;
                (e, rest) = rest.split_first_mut().unwrap();
            }
            if the_max == max {
                self.0.pop();
            } else {
                let last = self.0.last_mut().unwrap();
                last.start = last.end;
                last.end = max;
            }
            if the_min > min {
                self.0.insert(0, min..the_min);
            }
        }
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    pub fn iter(&'a self) -> std::slice::Iter<'a, Range<T>> {
        self.0.iter()
    }
}

#[test]
fn test_sort() {
    let mut i = Intervals::from_tuples(&[(1u32, 4), (12, 13), (5, 10)]);
    let j = Intervals::from_tuples(&[(1u32, 4), (5, 10), (12, 13)]);
    i.sort();
    assert!(i.is_equal(&j));
}

#[test]
fn test_merge() {
    let mut i1 = Intervals::from_tuples(&[(1u32, 3), (3, 4), (5, 10), (12, 13)]);
    let mut i2 = Intervals::from_tuples(&[(1u32, 3), (1, 4), (3, 4), (5, 10), (12, 13)]);
    let mut i3 = Intervals::from_tuples(&[(1u32, 3), (1, 4), (3, 4), (5, 6), (5, 10), (12, 13)]);
    let j = Intervals::from_tuples(&[(1u32, 4), (5, 10), (12, 13)]);

    i1.merge();
    i2.merge();
    i3.merge();

    assert!(i1.is_equal(&j));
    assert!(i2.is_equal(&j));
    assert!(i3.is_equal(&j));
}

#[test]
fn test_complement() {
    let mut i = Intervals::from_tuples(&[(1u32, 4), (5, 10), (12, 13)]);
    let mut i2 = Intervals::from_tuples(&[(1u32, 3), (2, 4), (5, 10), (12, 13)]);
    let j1 = Intervals::from_tuples(&[(0u32, 1), (4, 5), (10, 12), (13, 20)]);

    i.complement(0, 20);
    i2.complement(0, 20);
    assert!(i.is_equal(&j1));
    assert!(i2.is_equal(&j1));
}
