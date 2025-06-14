//! A simple and generic implementation of an immutable interval tree.
//!
//! Original version is from
//! <https://github.com/main--/rust-intervaltree/blob/master/src/lib.rs>
//!
//! Modifications by Bing Guo:
//! 23-April-3: allow reusing interval vector (add method `new` and `clear_and_fill_with_iter`)

// #[cfg(not(feature = "std"))]
extern crate alloc;
// #[cfg(feature = "serde")]
extern crate serde;
// #[cfg(feature = "std")]
extern crate std;

// #[cfg(not(feature = "std"))]
use alloc::vec::{IntoIter, Vec};
use core::cmp;
use core::fmt::{Debug, Formatter, Result as FmtResult};
use core::iter::FromIterator;
use core::ops::Range;
use core::slice::Iter;
// #[cfg(feature = "serde")]
// use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
// #[cfg(feature = "std")]
// use std::vec::{IntoIter, Vec};

/// An element of an interval tree.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
// #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Element<K, V> {
    /// The range associated with this element.
    pub range: Range<K>,
    /// The value associated with this element.
    pub value: V,
}

impl<K, V> From<(Range<K>, V)> for Element<K, V> {
    fn from(tup: (Range<K>, V)) -> Element<K, V> {
        let (range, value) = tup;
        Element { range, value }
    }
}

#[derive(Clone, Debug, Hash)]
// #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Node<K, V> {
    pub element: Element<K, V>,
    pub max: K,
}

/// A simple and generic implementation of an immutable interval tree.
///
/// To build it, always use `FromIterator`. This is not very optimized
/// as it takes `O(log n)` stack (it uses recursion) but runs in `O(n log n)`.
#[derive(Clone, Debug, Hash)]
// #[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct IntervalTree<K, V> {
    data: Vec<Node<K, V>>,
}

impl<K: Ord + Clone, V> IntervalTree<K, V> {
    /// Create an empty object
    ///
    pub fn new(capacity: usize) -> Self {
        Self {
            data: Vec::<Node<K, V>>::with_capacity(capacity),
        }
    }

    /// build trees from nodes. this is more flexible than
    /// building from iterators.
    ///
    /// Note: the max value in the node is not used. It will be recalculated
    pub fn new_from_nodes(mut nodes: Vec<Node<K, V>>) -> Self {
        nodes
            .iter_mut()
            .for_each(|n| n.max = n.element.range.end.clone());
        nodes.sort_unstable_by(|a, b| a.element.range.start.cmp(&b.element.range.start));
        if !nodes.is_empty() {
            Self::update_max(&mut nodes);
        }

        IntervalTree { data: nodes }
    }

    /// extract nodes (buffer) to reduce allocation for building new trees
    ///
    /// The extracted nodes can be used from `new_from_nodes`
    pub fn into_nodes(&mut self) -> Vec<Node<K, V>> {
        std::mem::take(&mut self.data)
    }

    /// Clear interval vector and fill with value from iterator
    ///
    pub fn clear_and_fill_with_iter<T: IntoIterator<Item = I>, I: Into<Element<K, V>>>(
        &mut self,
        iter: T,
    ) {
        self.data.clear();
        self.data
            .extend(iter.into_iter().map(|i| i.into()).map(|element| Node {
                max: element.range.end.clone(),
                element,
            }));

        self.data
            .sort_unstable_by(|a, b| a.element.range.start.cmp(&b.element.range.start));

        if !self.data.is_empty() {
            Self::update_max(&mut self.data);
        }
    }
}

impl<K: Ord + Clone, V, I: Into<Element<K, V>>> FromIterator<I> for IntervalTree<K, V> {
    fn from_iter<T: IntoIterator<Item = I>>(iter: T) -> Self {
        let mut nodes: Vec<_> = iter
            .into_iter()
            .map(|i| i.into())
            .map(|element| Node {
                max: element.range.end.clone(),
                element,
            })
            .collect();

        nodes.sort_unstable_by(|a, b| a.element.range.start.cmp(&b.element.range.start));

        if !nodes.is_empty() {
            Self::update_max(&mut nodes);
        }

        IntervalTree { data: nodes }
    }
}

/// An iterator over all the elements in the tree (in no particular order).
pub struct TreeIter<'a, K: 'a, V: 'a>(Iter<'a, Node<K, V>>);

impl<'a, K: 'a, V: 'a> Iterator for TreeIter<'a, K, V> {
    type Item = &'a Element<K, V>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|x| &x.element)
    }
}

impl<'a, K: 'a + Ord, V: 'a> IntoIterator for &'a IntervalTree<K, V> {
    type Item = &'a Element<K, V>;
    type IntoIter = TreeIter<'a, K, V>;

    fn into_iter(self) -> TreeIter<'a, K, V> {
        self.iter()
    }
}

/// An iterator that moves out of an interval tree.
pub struct TreeIntoIter<K, V>(IntoIter<Node<K, V>>);

impl<K, V> IntoIterator for IntervalTree<K, V> {
    type Item = Element<K, V>;
    type IntoIter = TreeIntoIter<K, V>;

    fn into_iter(self) -> TreeIntoIter<K, V> {
        TreeIntoIter(self.data.into_iter())
    }
}

impl<K, V> Iterator for TreeIntoIter<K, V> {
    type Item = Element<K, V>;

    fn next(&mut self) -> Option<Element<K, V>> {
        self.0.next().map(|x| x.element)
    }
}

impl<K: Ord + Clone, V> IntervalTree<K, V> {
    fn update_max(nodes: &mut [Node<K, V>]) -> K {
        assert!(!nodes.is_empty());
        let i = nodes.len() / 2;
        if nodes.len() > 1 {
            {
                let (left, rest) = nodes.split_at_mut(i);
                if !left.is_empty() {
                    rest[0].max = cmp::max(rest[0].max.clone(), Self::update_max(left));
                }
            }

            {
                let (rest, right) = nodes.split_at_mut(i + 1);
                if !right.is_empty() {
                    rest[i].max = cmp::max(rest[i].max.clone(), Self::update_max(right));
                }
            }
        }

        nodes[i].max.clone()
    }
}

impl<K: Ord, V> IntervalTree<K, V> {
    fn todo(&self) -> TodoVec {
        let mut todo = SmallVec::new();
        if !self.data.is_empty() {
            todo.push((0, self.data.len()));
        }
        todo
    }

    /// Queries the interval tree for all elements overlapping a given interval.
    ///
    /// This runs in `O(log n + m)`.
    pub fn query(&self, range: Range<K>) -> QueryIter<K, V> {
        QueryIter {
            todo: self.todo(),
            tree: self,
            query: Query::Range(range),
        }
    }

    /// Queries the interval tree for all elements containing a given point.
    ///
    /// This runs in `O(log n + m)`.
    pub fn query_point(&self, point: K) -> QueryIter<K, V> {
        QueryIter {
            todo: self.todo(),
            tree: self,
            query: Query::Point(point),
        }
    }

    /// Returns an iterator over all elements in the tree (in no particular order).
    pub fn iter(&self) -> TreeIter<K, V> {
        TreeIter(self.data.iter())
    }

    /// Returns an iterator over all elements in the tree, sorted by `Element.range.start`.
    ///
    /// This is currently identical to `IntervalTree::iter` because the internal structure
    /// is already sorted this way, but may not be in the future.
    pub fn iter_sorted(&self) -> impl Iterator<Item = &Element<K, V>> {
        TreeIter(self.data.iter())
    }
}

#[derive(Clone)]
enum Query<K> {
    Point(K),
    Range(Range<K>),
}

impl<K: Ord> Query<K> {
    fn point(&self) -> &K {
        match *self {
            Query::Point(ref k) => k,
            Query::Range(ref r) => &r.start,
        }
    }

    fn go_right(&self, start: &K) -> bool {
        match *self {
            Query::Point(ref k) => k >= start,
            Query::Range(ref r) => &r.end > start,
        }
    }

    fn intersect(&self, range: &Range<K>) -> bool {
        match *self {
            Query::Point(ref k) => k < &range.end,
            Query::Range(ref r) => r.end > range.start && r.start < range.end,
        }
    }
}

type TodoVec = SmallVec<[(usize, usize); 16]>;

/// Iterator for query results.
pub struct QueryIter<'a, K: 'a, V: 'a> {
    tree: &'a IntervalTree<K, V>,
    todo: TodoVec,
    query: Query<K>,
}

impl<K: Ord + Clone, V> Clone for QueryIter<'_, K, V> {
    fn clone(&self) -> Self {
        QueryIter {
            tree: self.tree,
            todo: self.todo.clone(),
            query: self.query.clone(),
        }
    }
}

impl<K: Ord + Clone + Debug, V: Debug> Debug for QueryIter<'_, K, V> {
    fn fmt(&self, fmt: &mut Formatter) -> FmtResult {
        let v: Vec<_> = (*self).clone().collect();
        write!(fmt, "{:?}", v)
    }
}

impl<'a, K: Ord, V> Iterator for QueryIter<'a, K, V> {
    type Item = &'a Element<K, V>;

    fn next(&mut self) -> Option<&'a Element<K, V>> {
        while let Some((s, l)) = self.todo.pop() {
            let i = s + l / 2;

            let node = &self.tree.data[i];
            if self.query.point() < &node.max {
                // push left
                {
                    let leftsz = i - s;
                    if leftsz > 0 {
                        self.todo.push((s, leftsz));
                    }
                }

                if self.query.go_right(&node.element.range.start) {
                    // push right
                    {
                        let rightsz = l + s - i - 1;
                        if rightsz > 0 {
                            self.todo.push((i + 1, rightsz));
                        }
                    }

                    // finally, search this
                    if self.query.intersect(&node.element.range) {
                        return Some(&node.element);
                    }
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::iter;

    fn verify(tree: &IntervalTree<u32, u32>, i: u32, expected: &[u32]) {
        let mut v1: Vec<_> = tree.query_point(i).map(|x| x.value).collect();
        v1.sort();
        let mut v2: Vec<_> = tree.query(i..(i + 1)).map(|x| x.value).collect();
        v2.sort();
        assert_eq!(v1, expected);
        assert_eq!(v2, expected);
    }

    #[test]
    fn it_works() {
        let tree: IntervalTree<u32, u32> = [
            (0..3, 1),
            (1..4, 2),
            (2..5, 3),
            (3..6, 4),
            (4..7, 5),
            (5..8, 6),
            (4..5, 7),
            (2..7, 8),
        ]
        .iter()
        .cloned()
        .collect();

        verify(&tree, 0, &[1]);
        verify(&tree, 1, &[1, 2]);
        verify(&tree, 2, &[1, 2, 3, 8]);
        verify(&tree, 3, &[2, 3, 4, 8]);
        verify(&tree, 4, &[3, 4, 5, 7, 8]);
        verify(&tree, 5, &[4, 5, 6, 8]);
        verify(&tree, 6, &[5, 6, 8]);
        verify(&tree, 7, &[6]);
        verify(&tree, 8, &[]);
        verify(&tree, 9, &[]);
    }

    #[test]
    fn empty() {
        let tree: IntervalTree<u32, u32> = iter::empty::<Element<u32, u32>>().collect();
        verify(&tree, 42, &[]);
    }
}
