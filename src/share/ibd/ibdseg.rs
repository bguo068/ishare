use crate::gmap::GeneticMap;

/// struct representing an *encoded* IBD segment:
/// - for i and j, the lower 2 bits encode haplotype id, upper bits encodes individual id
/// - for s and e, the values are either chromosome-level bp poistion or genome-wide bp position
#[derive(Debug, PartialEq, PartialOrd, Copy, Clone)]
pub struct IbdSeg {
    pub i: u32,
    pub j: u32,
    pub s: u32,
    pub e: u32,
}

impl IbdSeg {
    /// Create a new IBD segment:
    /// - `i`: index of individual 1
    /// - `j`: index of individual 2
    /// - `m`: indicating which happlotype of individual 1
    /// - `n`: indicating which happlotype of individual 2
    /// - `s`: start position in base pair of the IBD segment (0-based)
    /// - `e`: end position in base pair of the IBD segment (0-based)
    /// - `pos_shift`: is an offset value that will be added `s` and `e`.
    ///
    /// Note 1:
    ///
    /// IBD segment start and end poisitionscan hold information of chromosome index by
    /// representing genome-wide coordinates (assuming all chromsomes are concatenated
    /// into a single one). See `GenomeInfo` methods
    /// [crate::genome::GenomeInfo::to_chr_pos] and
    /// [crate::genome::GenomeInfo::to_gw_pos] and
    /// for chromosomal and genome-wide position conversion.
    ///
    /// Note 2:
    ///
    /// The pair of haplotype indices can be 4 combination for original IBD segments
    /// - (1, 1)
    /// - (1, 2)
    /// - (2, 1)
    /// - (2, 2)
    ///
    /// For segment from merging, there is only one possibility:
    /// - (0, 0)
    ///
    pub fn new(i: u32, m: u8, j: u32, n: u8, s: u32, e: u32, pos_shift: u32) -> Self {
        IbdSeg {
            i: i * 4 + m as u32,
            j: j * 4 + n as u32,
            s: s + pos_shift,
            e: e + pos_shift,
        }
    }

    /// Normalize the IBD segment record by swapping (`i` and `j`) if `i` > `j`.
    /// This is useful to uniquely identify haplotype/individual pairs.
    pub fn normalized(&mut self) {
        if self.i > self.j {
            std::mem::swap(&mut self.i, &mut self.j);
        }
    }

    /// return raw `i` and `j`
    pub fn haplotype_pair_int(&self) -> (u32, u32) {
        (self.i, self.j)
    }

    /// return unpacked haplotype pairs(ind1, hap1, ind2, hap2)
    pub fn haplotype_pair(&self) -> (u32, u8, u32, u8) {
        (
            self.i >> 2,
            self.i as u8 & 0x3,
            self.j >> 2,
            self.j as u8 & 0x3,
        )
    }

    /// return unpacked individual pairs(ind1, ind2)
    pub fn individual_pair(&self) -> (u32, u32) {
        (self.i >> 2, self.j >> 2)
    }

    /// return start and end postion (genomw-wide coordinates)
    pub fn coords(&self) -> (u32, u32) {
        (self.s, self.e)
    }

    /// return IBD segment length in centimorgans
    pub fn get_seg_len_cm(&self, gmap: &GeneticMap) -> f32 {
        gmap.get_cm(self.e) - gmap.get_cm(self.s)
    }

    /// check if IBD segment is resulting from merging of more one original segments
    pub fn is_from_merge(&self) -> bool {
        ((self.i & 0x3) == 0) && ((self.j & 0x3) == 0)
    }
    /// check validity of IBD segments
    /// - the combinaton of hap1 and hap2 values are valid
    /// - start position is smaller then end position
    pub fn is_valid(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;

        (m <= 2) && (n <= 2) && ((m == 0) == (n == 0)) && (self.s < self.e)
    }
}
