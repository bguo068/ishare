use crate::gmap::GeneticMap;
use serde::{Deserialize, Serialize};
use std::path::Path;
/// struct representing an *encoded* IBD segment:
/// - for i and j, the lower 2 bits encode haplotype id, upper bits encodes individual id
/// - for s and e, the values are either chromosome-level bp poistion or genome-wide bp position
#[derive(Debug, PartialEq, PartialOrd, Eq, Ord, Copy, Clone, Serialize, Deserialize)]
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
    /// - `m`: indicating which happlotype of individual 1 (can be 0 and 1)
    /// - `n`: indicating which happlotype of individual 2 (can be 0 and 1)
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
    /// - (0, 0)
    /// - (0, 1)
    /// - (1, 0)
    /// - (1, 1)
    ///
    /// For segment from merging, there is only one possibility:
    /// - (2, 2)
    ///
    /// For segment from haploid individuals, haplotype indcies should be 3
    /// for all individuals and all segments
    /// - (3, 3)
    pub fn new(i: u32, m: u8, j: u32, n: u8, s: u32, e: u32, pos_shift: u32) -> Self {
        IbdSeg {
            i: i * 4 + m as u32,
            j: j * 4 + n as u32,
            s: s + pos_shift,
            e: e + pos_shift,
        }
    }

    /// Normalize the IBD segment record by swapping (`i` and `j`) if `i` < `j`.
    /// This is useful to uniquely identify haplotype/individual pairs.
    pub fn normalized(&mut self) {
        if self.i < self.j {
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

    pub fn is_haploid_ibd(&self) -> bool {
        ((self.i & 0x3) == 3) && ((self.j & 0x3) == 3)
    }
    pub fn is_dipoid_ibd(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;
        (m <= 2) && (n <= 2) && ((m == 2) == (n == 2))
    }

    /// check if IBD segment is resulting from merging of more one original segments
    pub fn is_from_merge(&self) -> bool {
        ((self.i & 0x3) == 2) && ((self.j & 0x3) == 2)
    }
    /// check validity of IBD segments
    /// - the combinaton of hap1 and hap2 values are valid
    /// - start position is smaller then end position
    pub fn is_valid(&self) -> bool {
        let m = self.i & 0x3;
        let n = self.j & 0x3;
        let valid_diploid = (m <= 2) && (n <= 2) && ((m == 2) == (n == 2));
        let valid_haploid = (m == 3) && (n == 3);
        let valid_coords = self.s < self.e;

        valid_coords && (valid_diploid || valid_haploid)
    }
}

pub fn write_ibdseg_vec(v: &Vec<IbdSeg>, out: impl AsRef<Path>) {
    use std::io::Write;
    let mut file = std::fs::File::create(out.as_ref())
        .map(std::io::BufWriter::new)
        .expect(&format!(
            "cannot create file {}",
            out.as_ref().to_str().unwrap()
        ));
    let sz = v.len() as u64;
    file.write(&sz.to_le_bytes()).unwrap();

    for seg in v.iter() {
        file.write(&seg.i.to_le_bytes()).unwrap();
        file.write(&seg.j.to_le_bytes()).unwrap();
        file.write(&seg.s.to_le_bytes()).unwrap();
        file.write(&seg.e.to_le_bytes()).unwrap();
    }
}

pub fn read_ibdseg_vec(out: impl AsRef<Path>) -> Vec<IbdSeg> {
    let mut file = std::fs::File::open(out.as_ref())
        .map(std::io::BufReader::new)
        .expect(&format!(
            "cannot read file {}",
            out.as_ref().to_str().unwrap()
        ));
    let mut byte8 = [0u8; 8];
    let mut byte4 = [0u8; 4];
    use std::io::Read;
    file.read_exact(&mut byte8).unwrap();
    let sz = u64::from_le_bytes(byte8);
    let mut v = Vec::with_capacity(sz as usize);

    for _ in 0..sz {
        file.read_exact(&mut byte4).unwrap();
        let i = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).unwrap();
        let j = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).unwrap();
        let s = u32::from_le_bytes(byte4);
        file.read_exact(&mut byte4).unwrap();
        let e = u32::from_le_bytes(byte4);
        let seg = IbdSeg { i, j, s, e };
        v.push(seg);
    }
    v
}
