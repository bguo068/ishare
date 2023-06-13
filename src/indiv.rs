use arrow::array::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::arrow_writer::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::fmt::Debug;
use std::sync::Arc;

#[derive(Serialize, Deserialize, Debug)]
pub struct Individuals {
    vec: Vec<String>,
    map: HashMap<String, usize>,
}
impl Individuals {
    pub fn from_txt_file(p: impl AsRef<Path>) -> Individuals {
        let vec: Vec<_> = std::fs::File::open(p.as_ref())
            .map(BufReader::new)
            .unwrap()
            .lines()
            .map(|r| r.unwrap())
            .collect();

        let map: HashMap<String, usize> = vec
            .iter()
            .enumerate()
            .map(|(i, sample_ref)| (sample_ref.clone(), i))
            .collect();

        Self { vec, map }
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
