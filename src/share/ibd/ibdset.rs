use super::ibdseg::IbdSeg;
use crate::container::intervals::Intervals;
use crate::container::intervaltree::IntervalTree;
use crate::genome::GenomeInfo;
use crate::genotype::common::GenotypeMatrix;
use crate::gmap::GeneticMap;
use crate::indiv::{Individuals, PloidyConverter};
use crate::share::mat::NamedMatrix;
use crate::site::Sites;
use ahash::HashMap;
use itertools::Itertools;
use rayon::prelude::*;
use rust_htslib::bgzf;
use rust_htslib::tpool::ThreadPool;
use std::path::Path;
use IbdSetPloidyStatus::*;
use IbdSetSortStatus::*;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IbdSetPloidyStatus {
    Diploid,
    DiploidMerged,
    Haploid,
    UnknownPloidy,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IbdSetSortStatus {
    SortedByIndividaulPair,
    SortedByHaplotypePair,
    Unsorted,
}

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
    ploidy_status: IbdSetPloidyStatus,
    sort_status: IbdSetSortStatus,
}

impl<'a> IbdSet<'a> {
    /// Create an empty IBD set with meta infomation
    pub fn new(gmap: &'a GeneticMap, ginfo: &'a GenomeInfo, inds: &'a Individuals) -> Self {
        Self {
            ibd: vec![],
            gmap,
            ginfo,
            inds,
            ploidy_status: UnknownPloidy,
            sort_status: Unsorted,
        }
    }

    /// Create IBD set from parts
    pub fn update_parts(
        &mut self,
        ibd: Vec<IbdSeg>,
        ploidy_status: IbdSetPloidyStatus,
        sort_status: IbdSetSortStatus,
    ) {
        self.ibd = ibd;
        self.ploidy_status = ploidy_status;
        self.sort_status = sort_status;
    }

    pub fn into_parts(&mut self) -> (Vec<IbdSeg>, IbdSetPloidyStatus, IbdSetSortStatus) {
        let mut v = vec![];
        std::mem::swap(&mut self.ibd, &mut v);
        (v, self.ploidy_status, self.sort_status)
    }

    pub fn into_vec(self) -> Vec<IbdSeg> {
        self.ibd
    }

    pub fn infer_ploidy(&mut self) {
        if self.is_valid_haploid_ibd() {
            self.ploidy_status = Haploid;
        } else if self.is_valid_diploid_ibd() {
            if self.has_merged_ibd() {
                self.ploidy_status = DiploidMerged;
            } else {
                self.ploidy_status = Diploid;
            }
        } else {
            self.ploidy_status = UnknownPloidy;
        }
    }

    /// Faster way to check ploidy status
    ///
    /// Faster than using methods that iterate over the whole IBD set
    ///  - self.has_merged_ibd();
    ///  - self.is_valid_diploid_ibd();
    ///  - self.is_valid_haploid_ibd();
    ///
    pub fn get_ploidy_status(&self) -> IbdSetPloidyStatus {
        self.ploidy_status
    }

    pub fn infer_sort_status(&mut self) {
        if self.is_sorted_by_haplotypes() {
            self.sort_status = SortedByHaplotypePair;
        } else if self.is_sorted_by_samples() {
            self.sort_status = SortedByIndividaulPair;
        } else {
            self.sort_status = Unsorted;
        }
    }

    /// Faster way to check sorted status
    ///
    /// Faster then using methods that iterate over the whole IBD set
    ///  - self.is_sorted_by_haplotypes();
    ///  - self.is_sorted_by_samples();
    ///
    pub fn get_sort_status(&self) -> IbdSetSortStatus {
        self.sort_status
    }

    /// Add a segment into the set
    ///
    /// Note add segment put ibdset in unknow ploidy status and unsorted status
    /// Use infer_ploidy when IBD addition are completely done.
    pub fn add(&mut self, ibdseg: IbdSeg) {
        match self.ploidy_status {
            UnknownPloidy => {}
            _ => self.ploidy_status = UnknownPloidy,
        }
        match self.sort_status {
            Unsorted => {}
            _ => {
                self.sort_status = Unsorted;
            }
        }
        self.ibd.push(ibdseg)
    }

    /// return a reference of genetic map
    pub fn get_gmap(&self) -> &GeneticMap {
        self.gmap
    }
    /// return a reference of GenomeInfo
    pub fn get_ginfo(&self) -> &GenomeInfo {
        self.ginfo
    }
    /// return a reference of Individuals
    pub fn get_inds(&self) -> &Individuals {
        self.inds
    }

    /// read all `*.ibd.gz` file (in hap-IBD format) into the IBD set
    /// by globing the folder and calling [IbdSet::read_hapibd_file]
    pub fn read_hapibd_dir(&mut self, p: impl AsRef<Path>) {
        for entry in p.as_ref().read_dir().unwrap() {
            if let Ok(entry) = entry {
                let filename = entry.file_name();
                let filename = filename.as_os_str().to_str().unwrap();
                if !filename.ends_with("ibd.gz") {
                    continue;
                }
                let p = entry.path();
                self.read_hapibd_file(&p, None);
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
    pub fn read_hapibd_file(&mut self, p: impl AsRef<Path>, min_cm: Option<f32>) {
        use csv::ByteRecord;
        use csv::ReaderBuilder;

        let tpool = ThreadPool::new(10).unwrap();
        let mut reader = bgzf::Reader::from_path(p.as_ref()).unwrap();
        reader.set_thread_pool(&tpool).unwrap();

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(reader);

        let mut record = ByteRecord::new();
        let ind_map = self.inds.m();

        // let gw_chr_start_cm = self.gmap.get_gw_chr_start_cm_vec(self.ginfo);
        use std::str::from_utf8;

        println!("{}", p.as_ref().to_string_lossy());
        loop {
            // this allows skipping rows that has wrong number of fields
            match reader.read_byte_record(&mut record) {
                Ok(yes) => {
                    if !yes {
                        break;
                    }
                }
                e => {
                    eprintln!("Warning: {:?}", e);
                    continue;
                }
            }

            // filter out very short segments
            if let Some(min_cm) = min_cm {
                let cm: f32 = from_utf8(&record[7]).unwrap().parse::<f32>().unwrap();
                if cm < min_cm {
                    continue;
                }
            }
            // filter out segment is Ibd sample names is not in individuals
            let i = ind_map.get(from_utf8(&record[0]).unwrap());
            let j = ind_map.get(from_utf8(&record[2]).unwrap());
            if i.is_none() || j.is_none() {
                continue;
            }
            let i = *i.unwrap() as u32;
            let j = *j.unwrap() as u32;
            let m: u8 = match from_utf8(&record[1]).unwrap().parse::<u8>().unwrap() {
                1 => 0,
                2 => 1,
                0 => 2,
                _ => panic!(),
            };
            let n: u8 = match from_utf8(&record[3]).unwrap().parse::<u8>().unwrap() {
                1 => 0,
                2 => 1,
                0 => 2,
                _ => panic!(),
            };
            let chr = from_utf8(&record[4]).unwrap();
            let s = from_utf8(&record[5]).unwrap().parse::<u32>().unwrap() - 1;
            let e = from_utf8(&record[6]).unwrap().parse::<u32>().unwrap() - 1;

            let chrid = self.ginfo.idx[chr];
            let pos_shift = self.ginfo.gwstarts[chrid];

            let mut ibd = IbdSeg::new(i, m, j, n, s, e, pos_shift);
            ibd.normalized();
            assert!(ibd.is_valid());
            self.ibd.push(ibd);
            // println!("find ibd");
        }
    }

    /// read IBD segment from a file in `tskibd` format:
    /// - columns are delimited by `\t`
    /// - the first 7 columns are used
    ///     - Id1: individual 1
    ///     - Id2: individual 2
    ///     - Start
    ///     - End
    ///     - Ancestor
    ///     - Tmrca
    ///     - HasMutation
    /// - individual/chromomes are converted to index according to meta information
    /// - position are converted from 1-based to 0-based.
    /// - for this format, m/n are 3 to represent haploid genome
    pub fn read_tskibd_file(&mut self, p: impl AsRef<Path>, chrname: &str) {
        use csv::ReaderBuilder;
        use csv::StringRecord;

        let reader = bgzf::Reader::from_path(p.as_ref()).unwrap();

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_reader(reader);

        let mut record = StringRecord::new();
        let ind_map = self.inds.m();
        let chr_map = &self.get_ginfo().idx;
        let chrid = chr_map[chrname];
        let chrmsize = self.ginfo.chromsize[chrid];
        let pos_shift = self.ginfo.gwstarts[chrid];

        // let gw_chr_start_cm = self.gmap.get_gw_chr_start_cm_vec(self.ginfo);

        // check header column names
        let header = &reader.headers().unwrap();
        assert_eq!(&header[0], "Id1");
        assert_eq!(&header[1], "Id2");
        assert_eq!(&header[2], "Start");
        assert_eq!(&header[3], "End");
        assert_eq!(&header[4], "Ancestor");
        assert_eq!(&header[5], "Tmrca");
        assert_eq!(&header[6], "HasMutation");

        let mut counter_invalid_ibd = 0usize;
        while reader.read_record(&mut record).unwrap() {
            // haploid genome
            let m = 3;
            let n = 3;

            let i = match ind_map.get(&record[0]) {
                Some(id) => *id as u32,
                None => continue,
            };
            let j = match ind_map.get(&record[1]) {
                Some(id) => *id as u32,
                None => continue,
            };

            // allow converting floats to ints
            let mut s = record[2].parse::<f32>().unwrap() as u32;
            if s >= 1 {
                s -= 1;
            }
            let mut e = record[3].parse::<f32>().unwrap() as u32;
            assert!(e <= chrmsize, "e={e}, chrmsize: {chrmsize}");
            if e >= 1 {
                e -= 1;
            }
            let mut ibd = IbdSeg::new(i, m, j, n, s, e, pos_shift);
            ibd.normalized();
            if ibd.is_valid() {
                self.ibd.push(ibd);
            } else {
                counter_invalid_ibd += 1;
                println!("{:?}", ibd);
            }
        }
        if counter_invalid_ibd > 0 {
            eprintln!(
                "WARN: there are {counter_invalid_ibd} invalid records removed in file: {}",
                p.as_ref().to_str().unwrap()
            );
        }
    }

    /// read all `{chrname}.ibd` file (in hap-IBD format) into the IBD set
    /// by globing the folder and calling [IbdSet::read_hapibd_file]
    pub fn read_tskibd_dir(&mut self, p: impl AsRef<Path>) {
        for entry in p.as_ref().read_dir().unwrap() {
            if let Ok(entry) = entry {
                let filename = entry.file_name();
                let filename = filename.as_os_str().to_str().unwrap();
                if !filename.ends_with("ibd") {
                    continue;
                }
                let p = entry.path();
                let chrname = p.file_stem().unwrap().to_str().unwrap();
                self.read_tskibd_file(&p, chrname);
            }
        }
    }

    /// read all `*.ibd.gz` file (in hap-IBD format) into the IBD set
    /// by globing the folder and calling [IbdSet::read_hapibd_file]
    pub fn read_hmmibd_dir(&mut self, p: impl AsRef<Path>) {
        for entry in p.as_ref().read_dir().unwrap() {
            if let Ok(entry) = entry {
                let filename = entry.file_name();
                let filename = filename.as_os_str().to_str().unwrap();
                if !filename.ends_with(".hmm.txt") {
                    continue;
                }
                let p = entry.path();
                self.read_hmmibd_file(&p, Some(2.0));
            }
        }
    }

    pub fn read_hmmibd_file(&mut self, p: impl AsRef<Path>, min_seg_cm: Option<f32>) {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(p.as_ref())
            .unwrap();
        let mut rec = csv::ByteRecord::new();

        fn to_str(b: &[u8]) -> &str {
            std::str::from_utf8(b).unwrap()
        }

        let sam_map = self.inds.m();
        while reader.read_byte_record(&mut rec).unwrap() {
            //  sample1,
            //  sample2,
            //  chrname,
            //  start_pos,
            //  end_pos,
            //  ibd,
            //  n_snp,

            // skip nonibd information
            let ibd = to_str(&rec[5]).parse::<u8>().unwrap();
            if ibd != 0 {
                continue;
            }

            let chrname = to_str(&rec[2]);
            let chrid = self.ginfo.idx[chrname];
            let start_pos = to_str(&rec[3]).parse::<u32>().unwrap() - 1; // 1-based to 0-based
            let end_pos = to_str(&rec[4]).parse::<u32>().unwrap() - 1; // 1-based to 0-based
            let s = self.ginfo.to_gw_pos(chrid, start_pos);
            let e = s + end_pos - start_pos;

            // skip short ibd segments if a threshold is provided
            if let Some(min_cm) = min_seg_cm {
                let cm = self.gmap.get_cm_len(s, e);
                if cm < min_cm {
                    continue;
                }
            }

            // skip samples not in the individual map
            let i = match sam_map.get(to_str(&rec[0])) {
                Some(i) => *i as u32,
                None => continue,
            };
            let j = match sam_map.get(to_str(&rec[1])) {
                Some(j) => *j as u32,
                None => continue,
            };
            // set m/n =3 to indicate haploids.
            let m = 3;
            let n = 3;
            let pos_shift = 0; //already converted to genomew-wise position

            let mut ibdseg = IbdSeg::new(i, m, j, n, s, e, pos_shift);
            ibdseg.normalized();

            self.ibd.push(ibdseg);
        }
    }
    pub fn to_hapibd_file(&self) {
        todo!()
    }

    /// Sort IBD segment by individual pair indicies and then by coordinates
    pub fn sort_by_samples(&mut self) {
        self.ibd
            .par_sort_by_key(|s| (s.individual_pair(), s.coords()));
        self.sort_status = SortedByIndividaulPair;
    }

    pub fn is_sorted_by_samples(&self) -> bool {
        self.ibd
            .iter()
            .zip(self.ibd.iter().skip(1))
            .all(|(x, y)| (x.individual_pair(), x.coords()) <= (y.individual_pair(), y.coords()))
    }

    /// Sort IBD segment by haplotype pair indicies and then by coordinates
    pub fn sort_by_haplotypes(&mut self) {
        self.ibd
            .par_sort_by_key(|s| (s.haplotype_pair(), s.coords()));
        self.sort_status = SortedByHaplotypePair;
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
        let is_haploid = match self.ploidy_status {
            Haploid => true,
            _ => false,
        };
        assert!(is_haploid);

        self.inds = diploid_inds;
        self.ibd.iter_mut().for_each(|seg| {
            let (ind1, ind2) = seg.individual_pair();
            let (ind1, hap1) = ploidy_converter.h2d(ind1).unwrap();
            let (ind2, hap2) = ploidy_converter.h2d(ind2).unwrap();
            seg.i = (ind1 << 2) + hap1 as u32;
            seg.j = (ind2 << 2) + hap2 as u32;
            seg.normalized();
        });
        self.ploidy_status = Diploid;
    }
    pub fn covert_to_haploid(
        &mut self,
        haploid_inds: &'a Individuals,
        ploidy_converter: &PloidyConverter,
    ) {
        // check ploidy
        let is_diploid = match self.ploidy_status {
            Diploid => true,
            DiploidMerged => true,
            _ => false,
        };
        assert!(is_diploid);

        self.inds = haploid_inds;
        self.ibd.iter_mut().for_each(|seg| {
            let (ind1, hap1, ind2, hap2) = seg.haplotype_pair();
            let ind1 = ploidy_converter.d2h(ind1, hap1);
            let ind2 = ploidy_converter.d2h(ind2, hap2);
            let hap1 = 3;
            let hap2 = 3;
            seg.i = (ind1 << 2) + hap1 as u32;
            seg.j = (ind2 << 2) + hap2 as u32;
            seg.normalized();
        });

        self.ploidy_status = Haploid;
    }

    pub fn merge(&mut self) {
        let is_diploid_unmerged = match self.ploidy_status {
            Diploid => true,
            _ => false,
        };
        assert!(is_diploid_unmerged);

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
        self.ploidy_status = DiploidMerged;
    }

    pub fn merge_using_browning_method(
        &mut self,
        min_cm: f32,
        max_ndiscord: u32,
        gt: &GenotypeMatrix,
        site: &Sites,
    ) {
        let is_diploid_unmerged = match self.ploidy_status {
            Diploid => true,
            _ => false,
        };

        assert!(is_diploid_unmerged);
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
        self.ploidy_status = DiploidMerged;
    }

    pub fn iter(&self) -> impl Iterator<Item = &IbdSeg> {
        self.ibd.iter()
    }
    pub fn as_slice(&self) -> &[IbdSeg] {
        self.ibd.as_slice()
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
                // set value at symmetrical location
                mat.set_by_positions(ind2, ind1, tot);
            } else {
                let (ind1, hap1, ind2, hap2) = blk[0].haplotype_pair();
                assert!((hap1 == 0) || (hap1 == 1));
                assert!((hap2 == 0) || (hap2 == 1));
                let hid1 = ind1 * 2 + hap1 as u32;
                let hid2 = ind2 * 2 + hap2 as u32;
                let mut tot = 0.0f32;
                for r in blk {
                    let len = gmap.get_cm_len(r.s, r.e);
                    tot += len;
                }
                mat.set_by_positions(hid1, hid2, tot);
                // set value at symmetrical location
                mat.set_by_positions(hid2, hid1, tot);
            }
        }

        mat
    }

    pub fn is_valid_haploid_ibd(&self) -> bool {
        self.ibd.iter().all(|x| x.is_haploid_ibd())
    }
    pub fn is_valid_diploid_ibd(&self) -> bool {
        self.ibd.iter().all(|x| x.is_dipoid_ibd())
    }
    pub fn has_merged_ibd(&self) -> bool {
        self.ibd.iter().any(|x| x.is_from_merge())
    }

    pub fn filter_segments_by_min_cm(&mut self, min_cm: f64) {
        // mark short ibd segments for deletion
        for seg in self.ibd.as_mut_slice() {
            let cm = seg.get_seg_len_cm(self.gmap);
            if (cm as f64) < min_cm {
                seg.i = 0;
                seg.j = 0;
            }
        }
        // remove marked segments (retain unmarked segments)
        self.ibd.retain(|seg| (seg.i != 0) || (seg.j != 0))
    }

    /// Remove regions from each IBD segment
    ///
    /// if an IBD segment is contained with a region, the whole IBD segment
    /// will be removed; if the segment is not touch by any regions,
    /// the segment will be intact; if the segment is partly overlapped by any
    /// regions, the parts of segments outside those regions will be kept.
    ///
    /// The resulting IBD segment can be further filtered:
    /// when min_cm is not None, the result short IBD segment will NOT be
    /// added to the output ibdset.
    ///
    pub fn remove_regions(
        &self,
        regions: &Intervals<u32>,
        res_ibd: &mut IbdSet,
        min_cm: Option<f32>,
    ) {
        // get complement of regions
        let mut regions = regions.clone();
        let genome_size = self.get_ginfo().get_total_len_bp();
        regions.complement(0, genome_size);
        // generate ibdseg that intersect the complement regions
        let tree = IntervalTree::from_iter(regions.iter().map(|x| (x.to_owned(), ())));
        let gmap = self.gmap;
        for seg in self.iter() {
            for element in tree.query(seg.e..seg.i) {
                let mut seg = seg.clone();
                if seg.s < element.range.start {
                    seg.s = element.range.start;
                }
                if seg.e > element.range.end {
                    seg.e = element.range.end;
                }
                match min_cm {
                    // if provided a threshold and segment is too short,
                    // the segment will be no added to output
                    Some(min_cm) if seg.get_seg_len_cm(gmap) < min_cm => {}
                    _ => {
                        res_ibd.add(seg);
                    }
                }
            }
        }
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
