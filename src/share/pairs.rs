use std::hash::Hash;

use ahash::HashSet;

/// iterator chunks of pairs
///
/// This struct is to reduce memory allocation for vector to store all pairs
#[derive(Debug)]
pub struct PairChunkIter<T: Send + Copy + Hash + Eq + Ord + std::fmt::Debug> {
    /// genome id only exists in list 1
    v1: Vec<T>,
    /// genome id only exists in list 1
    v2: Vec<T>,
    /// genome id exists in both list 1 and list2
    v12: Vec<T>,
    /// step size
    sz: usize,
    idx: usize,
    nsteps: [usize; 3],
}

impl<T> PairChunkIter<T>
where
    T: Send + Copy + Hash + Eq + Ord + std::fmt::Debug,
{
    pub fn new(a: &[T], b: &[T], sz: usize) -> Self {
        let h1: HashSet<T> = a.iter().map(|t| *t).collect();
        let h2: HashSet<T> = b.iter().map(|t| *t).collect();

        let mut v1: Vec<_> = h1.difference(&h2).map(|x| *x).collect();
        v1.sort();
        let mut v2: Vec<_> = h2.difference(&h1).map(|x| *x).collect();
        v2.sort();
        let mut v12: Vec<_> = h2.intersection(&h1).map(|x| *x).collect();
        v12.sort();

        let calc_nstep = |i: usize| match i % sz {
            0 => i / sz,
            _ => i / sz + 1,
        };

        let nsteps = [
            calc_nstep(v1.len()),
            calc_nstep(v2.len()),
            calc_nstep(v12.len()),
        ];

        Self {
            v1,
            v2,
            v12,
            sz,
            idx: 0,
            nsteps,
        }
    }

    pub fn get_n_chunks(&self) -> usize {
        self.nsteps[0] * self.nsteps[1] + self.nsteps[2] * self.nsteps[2]
    }

    /// this ensure the pair chunk contains at most self.sz genomes
    ///
    /// Note: pairs include both (a, b) and (b, a) if a != b
    pub fn next_pair_chunks(
        &mut self,
        pairs: &mut Vec<(T, T)>,
        related: &mut Vec<T>,
        is_within: &mut bool,
    ) -> bool {
        pairs.clear();
        related.clear();
        if self.idx < self.nsteps[0] * self.nsteps[1] {
            let i = self.idx / self.nsteps[1];
            let j = self.idx % self.nsteps[1];
            for a in self.v1[(i * self.sz)..].iter().take(self.sz) {
                for b in self.v2[(j * self.sz)..].iter().take(self.sz) {
                    pairs.push((*a, *b));
                }
            }
            related.extend(self.v1[(i * self.sz)..].iter().take(self.sz));
            related.extend(self.v2[(j * self.sz)..].iter().take(self.sz));
            related.sort();
            self.idx += 1;
            *is_within = false;
            return true;
        }
        let mn = self.nsteps[0] * self.nsteps[1];
        if self.idx < mn + self.nsteps[2] * self.nsteps[2] {
            let i = (self.idx - mn) / self.nsteps[2];
            let j = (self.idx - mn) % self.nsteps[2];
            for a in self.v12[(i * self.sz)..].iter().take(self.sz) {
                for b in self.v12[(j * self.sz)..].iter().take(self.sz) {
                    pairs.push((*a, *b));
                }
            }
            if pairs.len() > 0 {
                if i == j {
                    related.extend(self.v12[(i * self.sz)..].iter().take(self.sz));
                } else {
                    related.extend(self.v12[(i * self.sz)..].iter().take(self.sz));
                    related.extend(self.v12[(j * self.sz)..].iter().take(self.sz));
                }
            }
            related.sort();
            self.idx += 1;
            *is_within = true;
            return true;
        }
        false
    }
}

#[test]
fn test_pair_chuk_iter() {
    let a = vec![0, 1, 2, 3, 4, 5];
    let b = vec![3, 4, 5, 6, 7, 8];
    let sz = 2;
    let mut pc_iter = PairChunkIter::new(&a, &b, sz);
    let mut pairs = vec![];
    let mut related = vec![];
    let mut is_within = false;
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(0, 6), (0, 7), (1, 6), (1, 7)]);
    assert_eq!(related, vec![0, 1, 6, 7]);
    assert_eq!(is_within, false);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(0, 8), (1, 8)]);
    assert_eq!(related, vec![0, 1, 8]);
    assert_eq!(is_within, false);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(2, 6), (2, 7)]);
    assert_eq!(related, vec![2, 6, 7]);
    assert_eq!(is_within, false);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(2, 8)]);
    assert_eq!(related, vec![2, 8]);
    assert_eq!(is_within, false);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(3, 3), (3, 4), (4, 3), (4, 4)]);
    assert_eq!(related, vec![3, 4]);
    assert_eq!(is_within, true);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(3, 5), (4, 5)]);
    assert_eq!(related, vec![3, 4, 5]);
    assert_eq!(is_within, true);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(5, 3), (5, 4)]);
    assert_eq!(related, vec![3, 4, 5]);
    assert_eq!(is_within, true);
    pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(pairs, vec![(5, 5)]);
    assert_eq!(is_within, true);
    assert_eq!(related, vec![5]);
    let ret = pc_iter.next_pair_chunks(&mut pairs, &mut related, &mut is_within);
    assert_eq!(ret, false);
}
