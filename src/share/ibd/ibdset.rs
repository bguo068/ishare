use super::ibdseg::IbdSeg;
use crate::container::intervals::Intervals;
use crate::genome::GenomeInfo;
use crate::genotype::common::GenotypeMatrix;
use crate::gmap::GeneticMap;
use crate::indiv::{Individuals, PloidyConverter};
use crate::share::mat::NamedMatrix;
use crate::site::Sites;
use itertools::Itertools;
use rust_htslib::bgzf;
use std::collections::HashMap;
use std::path::Path;

/// A struct representing a set of IBD segments with references to meta information
/// - GeneticMap
/// - GenomeInfo
/// - Individual info
///
pub struct IbdSet<'a> {
    ibd: Vec<IbdSeg>,
    gmap: &'a GeneticMap,
    ginfo: &'a GenomeInfo,
    inds: &'a Individuals,
}

impl<'a> IbdSet<'a> {
    /// Create an empty IBD set with meta infomation
    pub fn new(gmap: &'a GeneticMap, ginfo: &'a GenomeInfo, inds: &'a Individuals) -> Self {
        Self {
            ibd: vec![],
            gmap,
            ginfo,
            inds,
        }
    }

    /// Add a segment into the set
    pub fn add(&mut self, ibdseg: IbdSeg) {
        self.ibd.push(ibdseg)
    }

    /// return a reference of genetic map
    pub fn get_gmap(&self) -> &GeneticMap {
        self.gmap
    }

    /// read all `*.ibd.gz` file (in hap-IBD format) into the IBD set
    /// by globing the folder and calling [IbdSet::read_hapibd_file]
    pub fn read_hapibd_dir(&mut self, p: impl AsRef<Path>) {
        for entry in p.as_ref().read_dir().unwrap() {
            if let Ok(entry) = entry {
                if !entry.file_type().unwrap().is_file() {
                    continue;
                }
                let filename = entry.file_name();
                let filename = filename.as_os_str().to_str().unwrap();
                if !filename.ends_with("ibd.gz") {
                    continue;
                }
                let p = entry.path();
                self.read_hapibd_file(&p);
            }
        }
    }

    /// read IBD segment from a file in `hap-ibd` format:
    /// - columns are delimited by `\t`
    /// - the first 7 columns are used
    ///     - individual 1
    ///     - hapl
    ///     - individual 2
    ///     - hap2
    ///     - chromsome
    ///     - start
    ///     - end
    /// - individual/chromomes are converted to index according to meta information
    /// - position are converted from 1-based to 0-based.
    pub fn read_hapibd_file(&mut self, p: impl AsRef<Path>) {
        use csv::ReaderBuilder;
        use csv::StringRecord;

        let reader = bgzf::Reader::from_path(p.as_ref()).unwrap();

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);

        let mut record = StringRecord::new();
        let ind_map = self.inds.m();

        // let gw_chr_start_cm = self.gmap.get_gw_chr_start_cm_vec(self.ginfo);

        while reader.read_record(&mut record).unwrap() {
            let i = ind_map[&record[0]] as u32;
            let m: u8 = record[1].parse::<u8>().unwrap();
            let j = ind_map[&record[2]] as u32;
            let n: u8 = record[3].parse::<u8>().unwrap();
            let chr = &record[4];
            let s = record[5].parse::<u32>().unwrap() - 1;
            let e = record[6].parse::<u32>().unwrap() - 1;
            // let l: f32 = record[7].parse::<f32>().unwrap();

            let chrid = self.ginfo.idx[chr];
            let pos_shift = self.ginfo.gwstarts[chrid];

            let mut ibd = IbdSeg::new(i, m, j, n, s, e, pos_shift);
            ibd.normalized();
            assert!(ibd.is_valid());
            self.ibd.push(ibd);
        }
    }

    /// Sort IBD segment by individual pair indicies and then by coordinates
    pub fn sort_by_samples(&mut self) {
        self.ibd.sort_by_key(|s| (s.individual_pair(), s.coords()));
    }

    pub fn is_sorted_by_samples(&self) -> bool {
        self.ibd
            .iter()
            .zip(self.ibd.iter().skip(1))
            .all(|(x, y)| (x.individual_pair(), x.coords()) <= (y.individual_pair(), y.coords()))
    }

    /// Sort IBD segment by haplotype pair indicies and then by coordinates
    pub fn sort_by_haplotypes(&mut self) {
        self.ibd.sort_by_key(|s| (s.haplotype_pair(), s.coords()));
    }

    pub fn is_sorted_by_haplotypes(&self) -> bool {
        self.ibd
            .iter()
            .zip(self.ibd.iter().skip(1))
            .all(|(x, y)| (x.haplotype_pair(), x.coords()) <= (y.haplotype_pair(), y.coords()))
    }

    /// Check IbdSet 1 and 2 sharing the same individuals (also in the same order)
    /// This is to make sure individual indices between the two sets are
    /// comparatible.
    pub fn has_same_individuals(&self, other: &Self) -> bool {
        &self.inds.v() == &other.inds.v()
    }

    pub fn covert_to_het_diploid(
        &mut self,
        diploid_inds: &'a Individuals,
        ploidy_converter: &PloidyConverter,
    ) {
        self.inds = diploid_inds;
        self.ibd.iter_mut().for_each(|seg| {
            let (ind1, ind2) = seg.individual_pair();
            let (ind1, hap1) = ploidy_converter.h2d(ind1).unwrap();
            let (ind2, hap2) = ploidy_converter.h2d(ind2).unwrap();
            seg.i = (ind1 << 2) + hap1 as u32;
            seg.j = (ind2 << 2) + hap2 as u32;
            seg.normalized();
        });
    }
    pub fn covert_to_haploid(
        &mut self,
        haploid_inds: &'a Individuals,
        ploidy_converter: &PloidyConverter,
    ) {
        self.inds = haploid_inds;
        self.ibd.iter_mut().for_each(|seg| {
            let (ind1, hap1, ind2, hap2) = seg.haplotype_pair();
            let ind1 = ploidy_converter.d2h(ind1, hap1);
            let ind2 = ploidy_converter.d2h(ind2, hap2);
            let hap1 = 0;
            let hap2 = 0;
            seg.i = (ind1 << 2) + hap1 as u32;
            seg.j = (ind2 << 2) + hap2 as u32;
            seg.normalized();
        });
    }

    pub fn merge(&mut self) {
        self.sort_by_samples();
        self.ibd
            .iter_mut()
            .group_by(|x| x.individual_pair())
            .into_iter()
            .for_each(|(_pair, mut grp)| {
                let mut ibd1 = grp.next().unwrap();
                for (i, x) in grp.enumerate() {
                    if i == 0 {
                        ibd1 = x;
                    } else {
                        let ibd2 = x;
                        let (_s1, e1) = (ibd1.s, ibd1.e);
                        let (s2, e2) = (ibd2.s, ibd2.e);
                        // assert!(s1<=s2);

                        if e1 >= s2 {
                            // overlapping
                            if e2 > e1 {
                                ibd1.e = e2;
                            }
                            ibd2.i = 0;
                            ibd2.j = 0;
                        } else {
                            // nochange if eveything else
                            ibd1 = ibd2; // point to next record
                        }
                    }
                }
            });
        // clean up unused record
        self.ibd.retain(|seg| !((seg.i == 0) && (seg.j == 0)));
    }

    pub fn merge_using_browning_method(
        &mut self,
        min_cm: f32,
        max_ndiscord: u32,
        gt: &GenotypeMatrix,
        site: &Sites,
    ) {
        let ginfo = self.ginfo;
        let gmap = self.gmap;

        let pos_map: HashMap<u32, usize> = site
            .get_gw_pos_slice()
            .iter()
            .enumerate()
            .map(|(i, p)| (*p, i))
            .collect();

        self.sort_by_samples();
        self.ibd
            .iter_mut()
            .group_by(|x| x.individual_pair())
            .into_iter()
            .for_each(|(_pair, mut grp)| {
                let mut ibd1 = grp.next().unwrap();
                for (i, x) in grp.enumerate() {
                    if i == 0 {
                        ibd1 = x;
                    } else {
                        let ibd2 = x;
                        let (s1, e1) = (ibd1.s, ibd1.e);
                        let (s2, e2) = (ibd2.s, ibd2.e);

                        // assert!(s1<=s2);
                        let (chr1, _, _) = ginfo.to_chr_pos(s1);
                        let (chr2, _, _) = ginfo.to_chr_pos(s2);

                        if e1 >= s2 {
                            // overlapping
                            if e2 > e1 {
                                ibd1.e = e2;
                            }
                            ibd2.i = 0;
                            ibd2.j = 0;
                        } else if (chr1 == chr2)
                            && (gmap.get_cm(s2) - gmap.get_cm(e1) < min_cm)
                            && gt.has_too_many_discod_sites(
                                ibd1.i as usize,
                                ibd1.j as usize,
                                ibd2.i as usize,
                                ibd2.j as usize,
                                pos_map[&e1],
                                pos_map[&s2],
                                max_ndiscord,
                            )
                        {
                            // merge if overlapping or nonoverlapping and close and not too many discrod
                            //     update ibd1 and mark ibd2 as unused record
                            if e2 > e1 {
                                ibd1.e = e2;
                            }
                            ibd2.i = 0;
                            ibd2.j = 0;
                        } else {
                            // nochange if eveything else
                            ibd1 = ibd2; // point to next record
                        }
                    }
                }
            });
        // clean up unused record
        self.ibd.retain(|seg| !((seg.i == 0) && (seg.j == 0)));
    }

    pub fn iter(&self) -> impl Iterator<Item = &IbdSeg> {
        self.ibd.iter()
    }

    pub fn get_gw_total_ibd_matrix(&self, ignore_hap: bool) -> NamedMatrix<f32> {
        let gmap = self.gmap;
        let mut mat;
        if ignore_hap {
            assert!(self.is_sorted_by_samples());
            let nind = self.inds.v().len() as u32;
            mat = NamedMatrix::new_from_shape(nind, nind);
        } else {
            assert!(self.is_sorted_by_haplotypes());
            let nhap = self.inds.v().len() as u32 * 2;
            mat = NamedMatrix::new_from_shape(nhap, nhap);
        }

        let blockiter = IbdSetBlockIter::new(&self, ignore_hap);
        let mut itvs = Intervals::new();
        for blk in blockiter {
            if ignore_hap {
                itvs.clear();
                let it = blk.iter().map(|seg| (seg.s, seg.e));
                itvs.extend_from_iter(it);
                itvs.merge();
                let (ind1, ind2) = blk[0].individual_pair();
                let mut tot = 0.0f32;
                for r in itvs.iter() {
                    let len = gmap.get_cm_len(r.start, r.end);
                    tot += len;
                }
                mat.set_by_positions(ind1, ind2, tot);
            } else {
                let (ind1, hap1, ind2, hap2) = blk[0].haplotype_pair();
                assert!((hap1 == 1) || (hap1 == 2));
                assert!((hap2 == 1) || (hap2 == 2));
                let hap1 = hap1 as u32 - 1;
                let hap2 = hap2 as u32 - 1;
                let hid1 = ind1 * 2 + hap1;
                let hid2 = ind2 * 2 + hap2;
                let mut tot = 0.0f32;
                for r in blk {
                    let len = gmap.get_cm_len(r.s, r.e);
                    tot += len;
                }
                mat.set_by_positions(hid1, hid2, tot);
            }
        }

        mat
    }
}

pub struct IbdSetBlockIter<'a> {
    ibd: &'a [IbdSeg],
    ignore_hap: bool,
}

impl<'a> IbdSetBlockIter<'a> {
    pub fn new(ibd: &'a IbdSet, ignore_hap: bool) -> Self {
        Self {
            ibd: ibd.ibd.as_slice(),
            ignore_hap,
        }
    }

    pub fn peak(&self) -> Option<&IbdSeg> {
        self.ibd.first()
    }
}

impl<'a> Iterator for IbdSetBlockIter<'a> {
    type Item = &'a [IbdSeg];
    fn next(&mut self) -> Option<Self::Item> {
        if self.ibd.len() == 0 {
            None
        } else {
            let first = self.peak().unwrap();
            let e_opt = self.ibd.iter().position(|x| match self.ignore_hap {
                true => x.individual_pair() != first.individual_pair(),
                false => x.haplotype_pair_int() != first.haplotype_pair_int(),
            });
            let e = match e_opt {
                Some(end) => end,
                None => self.ibd.len(),
            };
            let (blk, rest) = self.ibd.split_at(e);
            self.ibd = rest;
            Some(blk)
        }
    }
}

pub struct IbdSetBlockPairIter<'a> {
    a: IbdSetBlockIter<'a>,
    b: IbdSetBlockIter<'a>,
    ignore_hap: bool,
}

impl<'a> IbdSetBlockPairIter<'a> {
    pub fn new(ibd1: &'a IbdSet, ibd2: &'a IbdSet, ignore_hap: bool) -> Self {
        let a = IbdSetBlockIter::new(ibd1, ignore_hap);
        let b = IbdSetBlockIter::new(ibd2, ignore_hap);
        Self { a, b, ignore_hap }
    }
}

impl<'a> Iterator for IbdSetBlockPairIter<'a> {
    type Item = (Option<&'a [IbdSeg]>, Option<&'a [IbdSeg]>);
    fn next(&mut self) -> Option<Self::Item> {
        let grp1 = self.a.peak().map(|x| match self.ignore_hap {
            false => x.haplotype_pair_int(),
            true => x.individual_pair(),
        });
        let grp2 = self.b.peak().map(|x| match self.ignore_hap {
            false => x.haplotype_pair_int(),
            true => x.individual_pair(),
        });
        match (grp1, grp2) {
            (Some(p1), Some(p2)) => {
                if p1 < p2 {
                    Some((self.a.next(), None))
                } else if p1 == p2 {
                    Some((self.a.next(), self.b.next()))
                } else {
                    Some((None, self.b.next()))
                }
            }
            (Some(_p1), None) => Some((self.a.next(), None)),
            (None, Some(_p2)) => Some((None, self.b.next())),
            _ => None,
        }
    }
}

pub struct IbdSetBlockIterMut<'a> {
    ibd: &'a mut [IbdSeg],
    ignore_hap: bool,
}

impl<'a> IbdSetBlockIterMut<'a> {
    pub fn new(ibd: &'a mut [IbdSeg], ignore_hap: bool) -> Self {
        Self { ibd, ignore_hap }
    }

    pub fn peek(&self) -> Option<&IbdSeg> {
        self.ibd.first()
    }
}

impl<'a> Iterator for IbdSetBlockIterMut<'a> {
    type Item = &'a mut [IbdSeg];
    fn next(&mut self) -> Option<Self::Item> {
        if self.ibd.len() == 0 {
            None
        } else {
            let first = self.peek().unwrap();
            let e_opt = self.ibd.iter().position(|x| match self.ignore_hap {
                true => x.individual_pair() != first.individual_pair(),
                false => x.haplotype_pair_int() != first.haplotype_pair_int(),
            });
            let e = match e_opt {
                Some(end) => end,
                None => self.ibd.len(),
            };
            let (blk, rest) = self.ibd.split_at_mut(e);

            let blk = unsafe { std::mem::transmute(blk) };
            let rest = unsafe { std::mem::transmute(rest) };
            self.ibd = rest;
            Some(blk)
        }
    }
}

pub struct IbdSetBlockPairIterMut<'a> {
    a: IbdSetBlockIterMut<'a>,
    b: IbdSetBlockIterMut<'a>,
    ignore_hap: bool,
}

impl<'a> IbdSetBlockPairIterMut<'a> {
    pub fn new(ibd1: &'a mut IbdSet, ibd2: &'a mut IbdSet, ignore_hap: bool) -> Self {
        let a = IbdSetBlockIterMut::new(&mut ibd1.ibd[..], ignore_hap);
        let b = IbdSetBlockIterMut::new(&mut ibd2.ibd[..], ignore_hap);
        Self { a, b, ignore_hap }
    }
}

impl<'a> Iterator for IbdSetBlockPairIterMut<'a> {
    type Item = (Option<&'a mut [IbdSeg]>, Option<&'a mut [IbdSeg]>);
    fn next(&mut self) -> Option<Self::Item> {
        let grp1 = self.a.peek().map(|x| match self.ignore_hap {
            false => x.haplotype_pair_int(),
            true => x.individual_pair(),
        });
        let grp2 = self.b.peek().map(|x| match self.ignore_hap {
            false => x.haplotype_pair_int(),
            true => x.individual_pair(),
        });
        match (grp1, grp2) {
            (Some(p1), Some(p2)) => {
                if p1 < p2 {
                    Some((self.a.next(), None))
                } else if p1 == p2 {
                    Some((self.a.next(), self.b.next()))
                } else {
                    Some((None, self.b.next()))
                }
            }
            (Some(_p1), None) => Some((self.a.next(), None)),
            (None, Some(_p2)) => Some((None, self.b.next())),
            _ => None,
        }
    }
}

#[test]
fn test_ibdset_iter() {
    let mut ginfo = GenomeInfo::new();
    ginfo.chromnames = vec!["chr1".to_owned(), "chr2".to_owned()];
    ginfo.chromsize = vec![100_000_000, 100_000_000];
    ginfo.gwstarts = vec![0, 100_000_000];
    ginfo
        .idx
        .extend(vec![("chr1".to_owned(), 0), ("chr2".to_owned(), 1)]);

    let gmap = GeneticMap::from_iter(vec![(0, 0.0), (200_000_000, 200.0)].into_iter());
    let inds = Individuals::from_iter(vec!["a", "b", "c", "d"].into_iter());

    let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
    let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);
    ibd1.add(IbdSeg::new(0, 0, 1, 0, 10, 100, 0));
    ibd1.add(IbdSeg::new(0, 0, 1, 0, 90, 110, 0));
    ibd1.add(IbdSeg::new(0, 0, 1, 1, 90, 110, 0));
    ibd1.add(IbdSeg::new(0, 0, 2, 0, 10, 100, 0));
    ibd1.add(IbdSeg::new(0, 0, 2, 0, 90, 110, 0));
    ibd1.add(IbdSeg::new(0, 0, 2, 1, 90, 110, 0));

    ibd2.add(IbdSeg::new(0, 0, 2, 0, 20, 100, 0));
    ibd2.add(IbdSeg::new(0, 0, 2, 0, 80, 110, 0));
    ibd2.add(IbdSeg::new(0, 0, 2, 1, 90, 130, 0));
    ibd2.add(IbdSeg::new(0, 0, 3, 0, 10, 100, 0));
    ibd2.add(IbdSeg::new(0, 0, 3, 0, 90, 110, 0));
    ibd2.add(IbdSeg::new(0, 0, 3, 1, 90, 110, 0));

    ibd1.sort_by_samples();
    ibd2.sort_by_samples();

    for (a, b) in IbdSetBlockPairIter::new(&ibd1, &ibd2, true) {
        match (a, b) {
            (Some(a), Some(b)) => {
                assert_eq!(a[0].individual_pair(), b[0].individual_pair());
            }
            _ => {}
        }
    }
    ibd1.sort_by_haplotypes();
    ibd2.sort_by_haplotypes();

    for (a, b) in IbdSetBlockPairIter::new(&ibd1, &ibd2, false) {
        match (a, b) {
            (Some(a), Some(b)) => {
                assert_eq!(a[0].haplotype_pair(), b[0].haplotype_pair());
            }
            _ => {}
        }
    }

    ibd1.sort_by_samples();
    ibd2.sort_by_samples();

    for (a, b) in IbdSetBlockPairIterMut::new(&mut ibd1, &mut ibd2, true) {
        match (a, b) {
            (Some(a), Some(b)) => {
                a[0].e += 1;
                assert_eq!(a[0].individual_pair(), b[0].individual_pair());
            }
            _ => {}
        }
    }
    ibd1.sort_by_haplotypes();
    ibd2.sort_by_haplotypes();

    for (a, b) in IbdSetBlockPairIterMut::new(&mut ibd1, &mut ibd2, false) {
        match (a, b) {
            (Some(a), Some(b)) => {
                a[0].e += 1;
                assert_eq!(a[0].haplotype_pair(), b[0].haplotype_pair());
            }
            _ => {}
        }
    }
}
