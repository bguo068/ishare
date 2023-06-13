use ishare::{
    container::{intervals::Intervals, intervaltree},
    genome, gmap,
    indiv::*,
    share::ibd::{ibdseg::*, ibdset::*},
};
use itertools;
use itertools::EitherOrBoth::*;
use itertools::Itertools;

fn main() {
    for set in
        "EPIGEN_BRAZIL_BAMBU,EPIGEN_BRAZIL_PELOTAS,EPIGEN_BRAZIL_SCAALA,full_PGP_LDGH".split(",")
    {
        println!("\n\n =================== {set} ========================\n");
        compare_ibd(set);
        compare_ibd_median_coverage_median_length(set);
    }
}

fn compare_ibd_median_coverage_median_length(set: &str) {
    // meta information
    let ginfo = genome::GenomeInfo::from_toml_file("genome.toml");
    let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
    let sample_fn = format!("testdata/{set}/samples.txt");
    println!("sample_fn = {sample_fn}");
    let inds = Individuals::from_txt_file(&sample_fn);

    // read / encode ibd segments
    let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
    let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);
    ibd1.read_hapibd_dir(&format!("testdata/{set}/gt/result/hapibd/"));
    ibd2.read_hapibd_dir(&format!("testdata/{set}/ph/result/hapibd/"));

    // sort by haplotype/sample pairs
    ibd1.sort_by_samples();
    ibd2.sort_by_samples();

    let gw_size_bp = ginfo.get_total_len_bp();

    let step = 100_000;

    let num_step = gw_size_bp / 100_000 + 1;
    let mut counter1 = vec![0usize; num_step as usize];
    let mut counter2 = vec![0usize; num_step as usize];

    for seg in ibd1.iter() {
        let s = ((seg.s - step / 2) / step + 1) as usize;
        let e = ((seg.e - step / 2) / step + 1) as usize;
        // println!("a={s}, b={e}");
        counter1[s..e].iter_mut().for_each(|x| *x += 1);
    }
    for seg in ibd2.iter() {
        let s = ((seg.s - step / 2) / step + 1) as usize;
        let e = ((seg.e - step / 2) / step + 1) as usize;
        counter2[s..e].iter_mut().for_each(|x| *x += 1);
    }

    counter1.sort();
    counter2.sort();

    let mut npairs = inds.v().len() * 2;
    npairs = npairs * (npairs - 1) / 2;
    let mut median1 = counter1[counter1.len() / 2] as f64;
    median1 = median1 / npairs as f64;

    let mut median2 = counter2[counter2.len() / 2] as f64;
    median2 = median2 / npairs as f64;

    println!("median ratio for GT: {}", median1);
    println!("median ratio for PH: {}", median2);
    println!("median ratio  PH/GT: {}", median2 / median1);

    let sum: f64 = counter1.iter().map(|x| *x as f64).sum();
    let n = counter1.len();
    let avg1 = sum / n as f64 / npairs as f64;
    let sum: f64 = counter2.iter().map(|x| *x as f64).sum();
    let n = counter2.len();
    let avg2 = sum / n as f64 / npairs as f64;

    println!("mean ratio for GT: {}", avg1);
    println!("mean ratio for PH: {}", avg2);
    println!("mean ratio  PH/GT: {}", avg2 / avg1);

    let mut length1: Vec<u32> = ibd1
        .iter()
        .map(|seg| (seg.get_seg_len_cm(&gmap) * 1000.0) as u32)
        .collect();
    length1.sort();

    let mut length2: Vec<u32> = ibd2
        .iter()
        .map(|seg| (seg.get_seg_len_cm(&gmap) * 1000.0) as u32)
        .collect();
    length2.sort();

    let median_len1 = length1[length1.len() / 2] as f64 / 1000.0;
    let median_len2 = length2[length2.len() / 2] as f64 / 1000.0;

    println!("median len cm for GT: {}", median_len1);
    println!("median len cm for PH: {}", median_len2);
    println!("median len cm  PH/GT: {}", median_len2 / median_len1);

    let sum1: f64 = ibd1
        .iter()
        .map(|seg| seg.get_seg_len_cm(&gmap) as f64)
        .sum();
    let avg1 = sum1 / ibd1.iter().count() as f64;
    let sum2: f64 = ibd2
        .iter()
        .map(|seg| seg.get_seg_len_cm(&gmap) as f64)
        .sum();
    let avg2 = sum2 / ibd2.iter().count() as f64;
    println!("mean len cm for GT: {}", avg1);
    println!("mean len cm for PH: {}", avg2);
    println!("mean len cm  PH/GT: {}", avg2 / avg1);
}

fn compare_ibd(set: &str) {
    // meta information
    let ginfo = genome::GenomeInfo::from_toml_file("genome.toml");
    let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
    let sample_fn = format!("testdata/{set}/samples.txt");
    let inds = Individuals::from_txt_file(&sample_fn);

    // read / encode ibd segments
    let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
    let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);
    ibd1.read_hapibd_dir(&format!("testdata/{set}/gt/result/hapibd/"));
    ibd2.read_hapibd_dir(&format!("testdata/{set}/ph/result/hapibd/"));

    // sort by haplotype/sample pairs
    ibd1.sort_by_samples();
    ibd2.sort_by_samples();

    // check if two ibd sets are comparable
    assert!(ibd1.has_same_individuals(&ibd2));

    // iter blocks of ibd records
    let it1 = ibd1
        .iter()
        .filter(|x| x.get_seg_len_cm(&gmap) >= 3.0)
        .group_by(|x| x.individual_pair());
    let it2 = ibd2
        .iter()
        .filter(|x| x.get_seg_len_cm(&gmap) >= 3.0)
        .group_by(|x| x.individual_pair());

    // iter paired blocks
    let it3 = itertools::merge_join_by(it1.into_iter(), it2.into_iter(), |a, b| a.0.cmp(&b.0));

    let mut left_only = 0u32;
    let mut right_only = 0u32;
    let mut both = 0u32;

    // some buffer variables
    let mut tree = intervaltree::IntervalTree::<u32, ()>::new(40);
    let mut x_intervals = Intervals::new();
    let mut y_intervals = Intervals::new();

    // counter/aggregate variables
    let mut a_overlapped_by_b_counts = [0usize; 7];
    let mut a_overlapped_by_b_sums = [0.0f32; 7];
    let mut b_overlapped_by_a_counts = [0usize; 7];
    let mut b_overlapped_by_a_sums = [0.0f32; 7];

    let mut pair_total = Vec::<(f32, f32)>::new();

    // compare blocks
    for x in it3 {
        match x {
            // ibd for a give sample pair only exists for IBD set A
            Left(x) => {
                left_only += 1;
                // flatten x
                flatten_into_intervals(&mut x_intervals, x.1);
                // println!("{:?}", x_intervals);

                // go over each query interval
                let mut totalibd = 0.0f32;
                for r in x_intervals.iter() {
                    let len = gmap.get_cm_len(r.start, r.end);
                    let c = cm_classifier(len);
                    a_overlapped_by_b_counts[c] += 1;
                    totalibd += len;
                }
                pair_total.push((totalibd, 0.0f32));
            }
            // ibd for a give sample pair exists for both IBD set A and B
            Both(x, y) => {
                // both
                // println!("hwy");
                both += 1;
                let mut totalibd_x = 0.0f32;
                let mut totalibd_y = 0.0f32;

                // flatten x
                flatten_into_intervals(&mut x_intervals, x.1);
                flatten_into_intervals(&mut y_intervals, y.1);
                // insert into tree
                tree.clear_and_fill_with_iter(y_intervals.iter().map(|x| (x.clone(), ())));

                // query tree
                for r in x_intervals.iter() {
                    let len = gmap.get_cm(r.end) - gmap.get_cm(r.start);
                    let mut overlap_len_total = 0.0f32;
                    for element in tree.query(r.clone()) {
                        let r2 = &element.range;
                        let mut s = r.start;
                        let mut e = r.end;
                        if s < r2.start {
                            s = r2.start;
                        }
                        if e > r2.end {
                            e = r2.end
                        }
                        let single_overlap_en = gmap.get_cm(e) - gmap.get_cm(s);
                        overlap_len_total += single_overlap_en;
                    }

                    assert!(len >= overlap_len_total);

                    let c = cm_classifier(len);
                    // println!("{c}, {len}");
                    a_overlapped_by_b_sums[c] += overlap_len_total / len; //ratio
                    a_overlapped_by_b_counts[c] += 1;
                    totalibd_x += len;
                }
                // insert into tree
                tree.clear_and_fill_with_iter(x_intervals.iter().map(|x| (x.clone(), ())));

                // query tree
                for r in y_intervals.iter() {
                    let len = gmap.get_cm(r.end) - gmap.get_cm(r.start);
                    let mut overlap_len_total = 0.0f32;
                    for element in tree.query(r.clone()) {
                        let r2 = &element.range;
                        let mut s = r.start;
                        let mut e = r.end;
                        if s < r2.start {
                            s = r2.start;
                        }
                        if e > r2.end {
                            e = r2.end
                        }
                        let single_overlap_en = gmap.get_cm(e) - gmap.get_cm(s);
                        overlap_len_total += single_overlap_en;
                    }

                    assert!(len >= overlap_len_total);
                    let c = cm_classifier(len);
                    b_overlapped_by_a_sums[c] += overlap_len_total / len;
                    b_overlapped_by_a_counts[c] += 1;
                    totalibd_y += len;
                }
                pair_total.push((totalibd_x, totalibd_y));
                // println!("Both: {:?} {:?}", x.0, y.0);
            }
            // ibd for a give sample pair only exists for IBD set B
            Right(y) => {
                // right
                right_only += 1;
                // println!("right: {:?}", y.0);
                // flatten x
                flatten_into_intervals(&mut y_intervals, y.1);

                // println!("{}", y_intervals.len());

                let mut totalibd = 0f32;
                for r in y_intervals.iter() {
                    let len = gmap.get_cm(r.end) - gmap.get_cm(r.start);
                    let c = cm_classifier(len);
                    b_overlapped_by_a_counts[c] += 1;
                    totalibd += len;
                }
                pair_total.push((0.0f32, totalibd));
            }
        }
    }

    println!("left_only={left_only}, right_only={right_only}, both={both}");
    println!(
        "all={}, predicted_all={}",
        left_only + both + right_only,
        (inds.v().len()) * (inds.v().len() - 1) / 2
    );

    println!(
        " -------------------------------------------
            A overlapped by B: 
            [2, 3):          {}
            [3, 4):          {}
            [4, 6):          {}
            [6, 10):         {}
            [10, 18):        {}
            [18, Inf):       {}
        ",
        a_overlapped_by_b_sums[0] / a_overlapped_by_b_counts[0] as f32,
        a_overlapped_by_b_sums[1] / a_overlapped_by_b_counts[1] as f32,
        a_overlapped_by_b_sums[2] / a_overlapped_by_b_counts[2] as f32,
        a_overlapped_by_b_sums[3] / a_overlapped_by_b_counts[3] as f32,
        a_overlapped_by_b_sums[4] / a_overlapped_by_b_counts[4] as f32,
        a_overlapped_by_b_sums[5] / a_overlapped_by_b_counts[5] as f32,
    );
    println!(
        " -------------------------------------------
            B overlapped by A: 
            [2, 3):          {}
            [3, 4):          {}
            [4, 6):          {}
            [6, 10):         {}
            [10, 18):        {}
            [18, Inf):       {}
        ",
        b_overlapped_by_a_sums[0] / b_overlapped_by_a_counts[0] as f32,
        b_overlapped_by_a_sums[1] / b_overlapped_by_a_counts[1] as f32,
        b_overlapped_by_a_sums[2] / b_overlapped_by_a_counts[2] as f32,
        b_overlapped_by_a_sums[3] / b_overlapped_by_a_counts[3] as f32,
        b_overlapped_by_a_sums[4] / b_overlapped_by_a_counts[4] as f32,
        b_overlapped_by_a_sums[5] / b_overlapped_by_a_counts[5] as f32,
    );

    ////////////////////////
    let mut counts = [0usize; 16];
    let mut total = 0usize;

    let gw_sz_cm = gmap.get_size_cm();
    for (t1, t2) in pair_total.iter() {
        let id = error_classifier((t1 - t2) / gw_sz_cm);
        counts[id] += 1;
        total += 1;
    }
    let boundaries = [
        -1.0f32, -0.8f32, -0.5f32, -0.27f32, -0.09f32, -0.03f32, -0.01f32, 0.0f32, 0.01f32,
        0.03f32, 0.09f32, 0.27f32, 0.5f32, 0.8f32, 1.0f32, 2.0f32,
    ];
    println!("\n------------ total ibd error / genome size");
    for i in 0..(boundaries.len() - 1) {
        let s = boundaries[i];
        let e = boundaries[i + 1];
        println!(
            "[{s:10.4}, {e:10.6})             {:10.6}",
            (counts[i] as f32) / (total as f32)
        );
    }
    let total_ibd_fn = format!("totalibd_{set}");
    write_pair_total(pair_total, &total_ibd_fn);
}

fn cm_classifier(cm: f32) -> usize {
    let boundary = [2.0f32, 3.0, 4.0, 6.0, 10.0, 18.0, f32::MAX];
    boundary.partition_point(|x| *x <= cm) - 1
}
fn error_classifier(cm: f32) -> usize {
    let boundary = [
        -1.0f32, -0.8f32, -0.5f32, -0.27f32, -0.09f32, -0.03f32, -0.01f32, 0.0f32, 0.01f32,
        0.03f32, 0.09f32, 0.27f32, 0.5f32, 0.8f32, 1.0f32, 2.0f32,
    ];
    boundary.partition_point(|x| *x <= cm) - 1
}
fn flatten_into_intervals<'a>(
    intervals: &mut Intervals<u32>,
    block: impl IntoIterator<Item = &'a IbdSeg>,
) {
    intervals.clear();
    block.into_iter().for_each(|e| intervals.push(e.s..e.e));
    intervals.sort();
    intervals.merge();
}

fn write_pair_total(pairtotal: Vec<(f32, f32)>, p: impl AsRef<std::path::Path>) {
    use arrow::array::{ArrayRef, Float32Array};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::arrow_writer::ArrowWriter;
    use parquet::file::properties::WriterProperties;
    use std::fs::File;
    use std::sync::Arc;

    let t1 = Float32Array::from_iter_values(pairtotal.iter().map(|x| x.0));
    let t2 = Float32Array::from_iter_values(pairtotal.iter().map(|x| x.1));
    let batch = RecordBatch::try_from_iter(vec![
        ("gt", Arc::new(t1) as ArrayRef),
        ("ph", Arc::new(t2) as ArrayRef),
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
