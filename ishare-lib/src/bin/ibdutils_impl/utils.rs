pub fn write_pair_total(
    pairtotal1: impl Iterator<Item = f32>,
    pairtotal2: impl Iterator<Item = f32>,
    colname1: &str,
    colname2: &str,
    p: impl AsRef<std::path::Path>,
) {
    use arrow::array::{ArrayRef, Float32Array};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::arrow_writer::ArrowWriter;
    use parquet::file::properties::WriterProperties;
    use std::fs::File;
    use std::sync::Arc;

    let t1 = Float32Array::from_iter_values(pairtotal1);
    let t2 = Float32Array::from_iter_values(pairtotal2);
    let batch = RecordBatch::try_from_iter(vec![
        (colname1.to_owned(), Arc::new(t1) as ArrayRef),
        (colname2.to_owned(), Arc::new(t2) as ArrayRef),
    ])
    .unwrap();

    let file = File::create(p.as_ref()).unwrap();

    // Default writer properties
    let props = WriterProperties::builder().build();

    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props)).unwrap();

    writer.write(&batch).expect("Writing batch");

    // writer must be closed to write footer
    writer.close().unwrap();
}
