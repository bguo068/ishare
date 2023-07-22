use crate::{
    container::intervaltree::IntervalTree,
    genome::GenomeInfo,
    gmap::GeneticMap,
    indiv::Individuals,
    share::ibd::{
        ibdseg::IbdSeg,
        ibdset::{IbdSet, IbdSetBlockIter},
    },
};

use super::fb::LASet;

#[derive(Clone, Copy, Debug)]
pub struct ASIbdSeg {
    pub ibd: IbdSeg,
    pub anc1: u8,
    pub anc2: u8,
}

pub struct ASIBDSet<'a> {
    asibd: Vec<ASIbdSeg>,
    gmap: &'a GeneticMap,
    ginfo: &'a GenomeInfo,
    inds: &'a Individuals,
    ancs: &'a [String],
}

impl<'a> ASIBDSet<'a> {
    /// Create an empty IBD set with meta infomation
    pub fn new(
        gmap: &'a GeneticMap,
        ginfo: &'a GenomeInfo,
        inds: &'a Individuals,
        ancs: &'a [String],
    ) -> Self {
        Self {
            asibd: vec![],
            gmap,
            ginfo,
            inds,
            ancs,
        }
    }
    pub fn add(&mut self, asibdseg: ASIbdSeg) {
        self.asibd.push(asibdseg)
    }

    pub fn get_asibd_from_ibdsets_and_laset(&mut self, ibdset: &IbdSet<'a>, la_set: &LASet) {
        let mut tree = IntervalTree::new(1000);
        for blk in IbdSetBlockIter::new(ibdset, false) {
            let (i, m, j, n) = blk[0].haplotype_pair();
            let hap1 = (i << 1) + m as u32;
            let hap2 = (j << 1) + n as u32;
            tree = la_set.get_hap_pair_la_segs2(hap1, hap2, tree);

            for ibdseg in blk {
                for elem in tree.query(ibdseg.s..ibdseg.e) {
                    let mut s = ibdseg.s;
                    let mut e = ibdseg.e;
                    if s < elem.range.start {
                        s = elem.range.start;
                    }
                    if e > elem.range.end {
                        e = elem.range.end;
                    }
                    let asibd = ASIbdSeg {
                        ibd: IbdSeg {
                            i: ibdseg.i,
                            j: ibdseg.j,
                            s,
                            e,
                        },
                        anc1: elem.value.0,
                        anc2: elem.value.1,
                    };
                    self.asibd.push(asibd);
                }
            }
        }
    }

    pub fn flush(&mut self, mut w: impl std::io::Write) {
        // use std::io::Write;
        for asibdseg in self.asibd.iter() {
            let ibd = &asibdseg.ibd;
            let (i, m, j, n) = ibd.haplotype_pair();
            let (_chrid, chrname, s) = self.ginfo.to_chr_pos(ibd.s);
            let e = s + (ibd.e - ibd.s);
            write!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                self.inds.v()[i as usize],
                m + 1,
                self.inds.v()[j as usize],
                n + 1,
                chrname,
                s + 1, // convert from 0-based to 1-based position
                e + 1, // convert from 0-based to 1-based position
                self.gmap.get_cm_len(ibd.s, ibd.e),
                self.ancs[asibdseg.anc1 as usize],
                self.ancs[asibdseg.anc2 as usize],
            )
            .unwrap();
        }
    }
}
