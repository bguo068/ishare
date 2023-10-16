use crate::container::intervaltree::{Element, IntervalTree};
use std::ops::Range;

pub struct CovCounter {
    tree: IntervalTree<u32, usize>,
    counts: Vec<usize>,
}
impl CovCounter {
    pub fn new(intervals: impl Iterator<Item = Range<u32>>) -> Self {
        let tree = IntervalTree::from_iter(
            intervals
                .enumerate()
                .map(|(i, r)| Element::from((r.clone(), i))),
        );
        let n = tree.iter().count();
        let counts = vec![0usize; n];
        CovCounter { tree, counts }
    }

    pub fn from_range(start: u32, end: u32, step: u32) -> Self {
        let it = (start..end)
            .step_by(step as usize)
            .enumerate()
            .map(|(i, x)| Element::from((x..(x + 1), i)));
        let tree = IntervalTree::from_iter(it);
        let counts = vec![0usize; tree.iter().count()];
        CovCounter { tree, counts }
    }
    pub fn count_over_interval(&mut self, r: &Range<u32>) {
        for e in self.tree.query(r.clone()) {
            self.counts[e.value] += 1;
        }
    }
    pub fn get_counts(&self) -> &[usize] {
        self.counts.as_slice()
    }

    pub fn get_intervals(&self) -> impl IntoIterator<Item = Range<u32>> + '_ {
        self.tree.iter_sorted().map(|x| x.range.start..x.range.end)
    }

    pub fn get_n_interval(&self) -> usize {
        self.tree.iter().count()
    }

    pub fn iter_sorted_start_end_count(&self) -> impl IntoIterator<Item = (u32, u32, usize)> + '_ {
        self.tree
            .iter_sorted()
            .map(|e| (e.range.start, e.range.end, self.counts[e.value]))
    }

    pub fn get_median_count(&self) -> f64 {
        let mut v = self.counts.clone();
        v.sort_unstable();

        let b = v[v.len() / 2] as f64;
        match v.len() % 2 == 0 {
            true => {
                let a = v[v.len() / 2 - 1] as f64;
                (a + b) / 2.0
            }
            false => b,
        }
    }

    pub fn get_mean_count(&self) -> f64 {
        let sum = self.counts.iter().sum::<usize>() as f64;
        let n = self.counts.len() as f64;
        sum / n
    }

    pub fn into_parquet(&mut self, p: impl AsRef<std::path::Path>) {
        use arrow::array::*;
        use arrow::record_batch::RecordBatch;
        use parquet::arrow::arrow_writer::ArrowWriter;
        use parquet::file::properties::WriterProperties;
        use std::sync::Arc;

        let starts: Vec<u32> = self.tree.iter_sorted().map(|x| x.range.start).collect();
        let starts = Arc::new(UInt32Array::from(starts)) as ArrayRef;

        let ends: Vec<u32> = self.tree.iter_sorted().map(|x| x.range.end).collect();
        let ends = Arc::new(UInt32Array::from(ends)) as ArrayRef;

        // when intervals used to construct CovCounter is not sorted
        // the counts need to be sorted by position
        let coverage: Vec<u64> = self
            .tree
            .iter_sorted()
            .map(|e| self.counts[e.value] as u64)
            .collect();

        let coverage = Arc::new(UInt64Array::from(coverage)) as ArrayRef;

        let batch = RecordBatch::try_from_iter(vec![
            ("Start", starts),
            ("End", ends),
            ("Coverage", coverage),
        ])
        .unwrap();

        let file = std::fs::File::create(p.as_ref()).unwrap();
        let prop = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(prop)).unwrap();

        writer.write(&batch).expect("Write arrow batch");
        writer.close().unwrap();
    }
}

#[test]
fn test_cov_counter() {
    let i = vec![10u32..11, 20..21, 30..31, 40..41];
    let mut counter = CovCounter::new(i.iter().map(|r| r.clone()));

    let j = vec![10u32..22, 10..21, 1..34];

    for r in j.iter() {
        counter.count_over_interval(r);
    }

    assert_eq!(counter.get_counts(), &[3, 3, 1, 0]);

    counter.into_parquet("tmp.cov.pq");
}
