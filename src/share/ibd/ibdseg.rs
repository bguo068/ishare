use crate::gmap::GeneticMap;

/// for i and j, the lower 2 bits encode haplotype id, upper bits;
/// for s and e, the values are either chromosome-level bp poistion or genome-wide bp position
#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
pub struct IbdSeg {
    pub i: u32,
    pub j: u32,
    pub s: u32,
    pub e: u32,
}

impl IbdSeg {
    pub fn new(i: u32, m: u8, j: u32, n: u8, s: u32, e: u32, pos_shift: u32) -> Self {
        IbdSeg {
            i: i * 4 + m as u32,
            j: j * 4 + n as u32,
            s: s + pos_shift,
            e: e + pos_shift,
        }
    }
    pub fn normalized(&mut self) {
        if self.i > self.j {
            std::mem::swap(&mut self.i, &mut self.j);
        }
    }

    pub fn haplotype_pair_int(&self) -> (u32, u32) {
        (self.i, self.j)
    }

    pub fn haplotype_pair(&self) -> (u32, u8, u32, u8) {
        (
            self.i >> 2,
            self.i as u8 & 0x3,
            self.j >> 2,
            self.j as u8 & 0x3,
        )
    }
    pub fn individual_pair(&self) -> (u32, u32) {
        (self.i >> 2, self.j >> 2)
    }

    pub fn coords(&self) -> (u32, u32) {
        (self.s, self.e)
    }
    pub fn get_seg_len_cm(&self, gmap: &GeneticMap) -> f32 {
        gmap.get_cm(self.e) - gmap.get_cm(self.s)
    }

    pub fn is_from_merge(&self) -> bool {
        ((self.i & 0x3) == 0) && ((self.j & 0x3) == 0)
    }
    pub fn is_valid(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;

        (m <= 2) && (n <= 2) && ((m == 0) == (n == 0)) && (self.s < self.e)
    }
}
