use super::ibdseg::IbdSeg;
use crate::genome::GenomeInfo;
use crate::genotype::common::GenotypeMatrix;
use crate::gmap::GeneticMap;
use crate::indiv::{Individuals, PloidyConverter};
use crate::site::Sites;
use itertools::Itertools;
use rust_htslib::bgzf;
use std::collections::HashMap;
use std::path::Path;

pub struct IbdSet<'a> {
    ibd: Vec<IbdSeg>,
    gmap: &'a GeneticMap,
    ginfo: &'a GenomeInfo,
    inds: &'a Individuals,
}

impl<'a> IbdSet<'a> {
    pub fn new(gmap: &'a GeneticMap, ginfo: &'a GenomeInfo, inds: &'a Individuals) -> Self {
        Self {
            ibd: vec![],
            gmap,
            ginfo,
            inds,
        }
    }
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

    pub fn sort_by_samples(&mut self) {
        self.ibd.sort_by_key(|s| (s.individual_pair(), s.coords()));
    }

    pub fn sort_by_haplotypes(&mut self) {
        self.ibd.sort_by_key(|s| (s.haplotype_pair(), s.coords()));
    }

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
}
