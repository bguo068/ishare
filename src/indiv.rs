use ahash::{HashMap, HashMapExt};
use arrow::array::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::path::Path;

use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::fmt::Debug;
use std::sync::Arc;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Individuals {
    vec: Vec<String>,
    map: HashMap<String, usize>,
}

pub enum PloidConvertDirection {
    Diploid2Haploid,
    Haploid2Diploid,
}

impl Individuals {
    /// Reads individual data from a text file.
    ///
    /// The function accepts three types of files:
    ///
    /// 1. Single Column:
    ///    Each line represents a name of a sample.
    /// 2. Two Columns (tab-separated):
    ///    This format indicates haploid to diploid conversion.  The first
    ///    column contains haploid sample names, and the second column contains
    ///    the corresponding diploid names.  Each diploid sample should have two
    ///    haploid samples mapped to it.  The first of the two haploid names
    ///    will be assigned haplotype index 1 (encoded as 0), and the second
    ///    haploid sample will be assigned haplotype index 2 (encoded as 1).
    /// 3. Three Columns (tab-separated):
    ///    This format indicates diploid to haploid conversion.  Each line
    ///    contains the diploid sample name followed by the names of the
    ///    converted haploid samples.  The second column contains the name of
    ///    haplotype 1 of the diploid sample, and the third column contains the
    ///    name of haplotype 2.
    ///
    /// Return Values:
    /// - The first `Individuals` is the "from" individuals.
    /// - `PloidConverter` contains two maps (`d2h` and `h2d`).
    /// - The second `Individuals` is the "to" individuals.
    /// - `PloidConvertDirection` indicates the conversion direction: either
    /// diploid-to-haploid or haploid-to-diploid.
    pub fn from_txt_file(
        p: impl AsRef<Path>,
    ) -> (
        Individuals,
        Option<(PloidyConverter, Individuals, PloidConvertDirection)>,
    ) {
        use std::io::read_to_string;
        let contents = read_to_string(File::open(p.as_ref()).unwrap()).unwrap();
        let ncolumns = contents.trim().lines().next().unwrap().split("\t").count();
        match ncolumns {
            1 => {
                let v: Vec<String> = contents.trim().lines().map(|x| x.to_owned()).collect();
                let m: HashMap<_, _> = v
                    .iter()
                    .enumerate()
                    .map(|(i, s)| (s.to_owned(), i))
                    .collect();
                (Self { vec: v, map: m }, None)
            }
            2 => {
                let direction = PloidConvertDirection::Haploid2Diploid;
                let mut v_h = vec![];
                let mut m_h = HashMap::new();
                let mut v_d = vec![];
                let mut m_d = HashMap::new();
                let mut h2dm = HashMap::new();
                let mut d2hm = HashMap::new();

                for line in contents.trim().lines() {
                    let mut fields = line.split("\t");
                    let h = fields.next().unwrap(); // column 1 is hap
                    let d = fields.next().unwrap(); // column 2 is dip
                    v_h.push(h.to_owned());
                    m_h.insert(h.to_owned(), m_h.len());

                    if !m_d.contains_key(d) {
                        // first hap
                        v_d.push(d.to_owned());
                        m_d.insert(d.to_owned(), m_d.len());
                        h2dm.insert(m_h[h] as u32, (m_d[d] as u32, 0));
                        d2hm.insert((m_d[d] as u32, 0), m_h[h] as u32);
                    } else {
                        // second hap

                        // assert hap1 is found ,hap2 not. This also prevents
                        // an errorneous third hap from being added.
                        assert!(d2hm.contains_key(&(m_d[d] as u32, 0)));
                        assert!(!d2hm.contains_key(&(m_d[d] as u32, 1)));

                        h2dm.insert(m_h[h] as u32, (m_d[d] as u32, 1));
                        d2hm.insert((m_d[d] as u32, 1), m_h[h] as u32);
                    }
                }

                (
                    Self { vec: v_h, map: m_h },
                    Some((
                        PloidyConverter { h2dm, d2hm },
                        Individuals { vec: v_d, map: m_d },
                        direction,
                    )),
                )
            }
            3 => {
                let direction = PloidConvertDirection::Diploid2Haploid;
                let mut v_d = vec![];
                let mut m_d = HashMap::new();
                let mut v_h = vec![];
                let mut m_h = HashMap::new();
                let mut h2dm = HashMap::new();
                let mut d2hm = HashMap::new();

                for line in contents.trim().lines() {
                    let mut fields = line.split("\t");
                    let d = fields.next().unwrap(); // column 1 is dip
                    let h1 = fields.next().unwrap(); // column 2 is hap
                    let h2 = fields.next().unwrap(); // column 2 is hap

                    v_d.push(d.to_owned());
                    m_d.insert(d.to_owned(), m_d.len());

                    v_h.push(h1.to_owned());
                    m_h.insert(h1.to_owned(), m_h.len());

                    v_h.push(h2.to_owned());
                    m_h.insert(h2.to_owned(), m_h.len());

                    h2dm.insert(m_h[h1] as u32, (m_d[d] as u32, 0));
                    d2hm.insert((m_d[d] as u32, 0), m_h[h1] as u32);

                    h2dm.insert(m_h[h2] as u32, (m_d[d] as u32, 1));
                    d2hm.insert((m_d[d] as u32, 1), m_h[h2] as u32);
                }

                (
                    Self { vec: v_d, map: m_d },
                    Some((
                        PloidyConverter { h2dm, d2hm },
                        Individuals { vec: v_h, map: m_h },
                        direction,
                    )),
                )
            }
            _ => panic!("samples/individuals file format error"),
        }
    }

    pub fn from_iter<'a>(it: impl Iterator<Item = &'a str> + 'a) -> Self {
        let mut v = Vec::<String>::new();
        let mut m = HashMap::<String, usize>::new();
        for e in it {
            m.insert(e.to_owned(), v.len());
            v.push(e.to_owned());
        }
        Self { vec: v, map: m }
    }
    pub fn v(&self) -> &Vec<String> {
        &self.vec
    }
    pub fn m(&self) -> &HashMap<String, usize> {
        &self.map
    }

    pub fn into_parquet_file(self, p: impl AsRef<Path>) {
        // array
        let v = StringArray::from(self.vec);
        let batch =
            RecordBatch::try_from_iter(vec![("chrnames", Arc::new(v) as ArrayRef)]).unwrap();
        // writer
        let file = File::create(p.as_ref()).unwrap();
        // -- default writer properties
        let props = WriterProperties::builder().build();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();
        // write batch
        writer.write(&batch).expect("Writing batch");
        // writer must be closed to write footer
        writer.close().unwrap();
    }

    pub fn from_parquet_file(p: impl AsRef<Path>) -> Self {
        let file = File::open(p).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
        let mut reader = builder.build().unwrap();

        let mut v = Vec::<String>::new();
        for record_batch in &mut reader {
            let record_batch = record_batch.unwrap();
            let it = record_batch
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .into_iter()
                .map(|x| x.unwrap().to_owned());

            v.extend(it);
        }

        let m: HashMap<_, _> = v
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_owned(), i))
            .collect();

        Individuals { vec: v, map: m }
    }

    pub fn get_ploidy_converter(&self) -> (Individuals, PloidyConverter) {
        let v: Vec<String> = self
            .v()
            .chunks_exact(2)
            .map(|s| format!("{}|{}", &s[0], &s[1]))
            .collect();
        let m = v.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect();
        let inds = Individuals { vec: v, map: m };

        // NOTE: if self.v().len() is an odd number, there will be data loss (leaving one sample out)
        let n = (self.v().len() / 2 * 2) as u32;
        let h2dm = (0..n).map(|x| (x, (x / 2, x as u8 % 2))).collect();
        let d2hm = (0..n).map(|x| ((x / 2, x as u8 % 2), x)).collect();
        let p = PloidyConverter { h2dm, d2hm };
        (inds, p)
    }
}

pub struct PloidyConverter {
    h2dm: HashMap<u32, (u32, u8)>,
    d2hm: HashMap<(u32, u8), u32>,
}

impl<'a> PloidyConverter {
    pub fn h2d(&self, h: u32) -> Option<(u32, u8)> {
        self.h2dm.get(&h).map(|(id, hid)| (*id, *hid))
    }
    pub fn d2h(&self, d: u32, which: u8) -> u32 {
        self.d2hm[&(d, which)]
    }
}
