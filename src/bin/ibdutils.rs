use clap::{Args, Parser, Subcommand};
use ishare::{
    genome::{self, GenomeInfo},
    gmap,
    indiv::*,
    share::ibd::{ibdseg::IbdSeg, ibdset::*, overlap},
};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(author, version, about, long_about=None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Encode genotype data from VCF file to tables or matrix (-m)
    Compare {
        /// Path to genome info toml file (input)
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        // Path to sample list file
        #[arg(short = 'S', long, default_value = "samples.txt")]
        samples: PathBuf,
        // Path to sample list file
        #[arg(short = 'o', long, default_value = "ibd_cmp_res.txt")]
        out: PathBuf,
        // fmt of ibd files, supported format 'hapibd', 'tskibd'
        #[arg(short = 'f', long, default_value = "hapibd")]
        fmt: String,
        // IBD directory 1
        ibd1_dir: PathBuf,
        // IBD directory 2
        ibd2_dir: PathBuf,
    },

    PlotIBD {
        #[arg(short = 'I', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        // Path to sample list file
        #[arg(short = 'S', long, default_value = "samples.txt")]
        samples: PathBuf,
        // Path to sample list file
        #[arg(short = 'o', long, default_value = "ibd_cmp_res.txt")]
        out: PathBuf,
        // fmt of ibd files, supported format 'hapibd', 'tskibd'
        #[arg(short = 'f', long, default_value = "hapibd")]
        fmt: String,
        // IBD directory 1
        ibd1_dir: PathBuf,
        // IBD directory 2
        ibd2_dir: PathBuf,

        #[command(flatten)]
        sample1: Sample1,
        #[command(flatten)]
        sample2: Sample2,
    },
}

#[derive(Args)]
#[group(required = true, multiple = false)]
struct Sample1 {
    #[arg(short = 'n', long)]
    ind_name1: Option<String>,
    #[arg(short = 'x', long)]
    ind_ix1: Option<u32>,
}

#[derive(Args)]
#[group(required = true, multiple = false)]
struct Sample2 {
    #[arg(short = 'N', long)]
    ind_name2: Option<String>,
    #[arg(short = 'X', long)]
    ind_ix2: Option<u32>,
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Compare {
            genome_info,
            samples,
            out,
            fmt,
            ibd1_dir,
            ibd2_dir,
        }) => {
            // files
            let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
            let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
            let inds = Individuals::from_txt_file(&samples);
            let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
            let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);

            if fmt == "hapibd" {
                ibd1.read_hapibd_dir(&ibd1_dir);
                ibd2.read_hapibd_dir(&ibd2_dir);
                ibd1.sort_by_samples();
                ibd2.sort_by_samples();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else if fmt == "tskibd" {
                ibd1.read_tskibd_dir(&ibd1_dir);
                ibd2.read_tskibd_dir(&ibd2_dir);
                ibd1.sort_by_haplotypes();
                ibd2.sort_by_haplotypes();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else if fmt == "hmmibd" {
                ibd1.read_hmmibd_dir(&ibd1_dir);
                ibd2.read_hmmibd_dir(&ibd2_dir);
                ibd1.sort_by_haplotypes();
                ibd2.sort_by_haplotypes();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else {
                panic!("format {} is not supported.", fmt);
            }

            {
                // overlapping analysis
                let oa = overlap::IbdOverlapAnalyzer::new(&mut ibd1, &mut ibd2, true);
                let res = oa.analzyze(Some(&[3.0f32, 4.0, 6.0, 10.0, 18.0]));
                res.to_csv(out);
            }
            {
                // total IBD analysis
                let mat1 = ibd1.get_gw_total_ibd_matrix(true);
                let mat2 = ibd2.get_gw_total_ibd_matrix(true);
                let n = inds.v().len();
                let pairs = (0..n - 1)
                    .map(|i| ((i + 1)..n).map(move |j| (i, j)))
                    .flatten();
                let it1 = pairs
                    .clone()
                    .map(|(i, j)| mat1.get_by_positions(i as u32, j as u32));
                let it2 = pairs
                    .clone()
                    .map(|(i, j)| mat2.get_by_positions(i as u32, j as u32));

                write_pair_total(it1, it2, out.with_extension("totibd"));
            }
        }
        Some(Commands::PlotIBD {
            genome_info,
            samples,
            out,
            fmt,
            ibd1_dir,
            ibd2_dir,
            sample1,
            sample2,
        }) => {
            let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
            let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
            let inds = Individuals::from_txt_file(&samples);
            let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
            let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);

            if fmt == "hapibd" {
                ibd1.read_hapibd_dir(&ibd1_dir);
                ibd2.read_hapibd_dir(&ibd2_dir);
                ibd1.sort_by_samples();
                ibd2.sort_by_samples();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else if fmt == "tskibd" {
                ibd1.read_tskibd_dir(&ibd1_dir);
                ibd2.read_tskibd_dir(&ibd2_dir);
                ibd1.sort_by_haplotypes();
                ibd2.sort_by_haplotypes();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else if fmt == "hmmibd" {
                ibd1.read_hmmibd_dir(&ibd1_dir);
                ibd2.read_hmmibd_dir(&ibd2_dir);
                ibd1.sort_by_haplotypes();
                ibd2.sort_by_haplotypes();
                ibd1.infer_ploidy();
                ibd2.infer_ploidy();
            } else {
                panic!("format {} is not supported.", fmt);
            }

            let mut id1 = match sample1.ind_ix1 {
                Some(id) => {
                    assert!(id < inds.v().len() as u32);
                    id
                }
                None => {
                    let s = sample1.ind_name1.as_ref().unwrap();
                    inds.m()[s] as u32
                }
            };
            let mut id2 = match sample2.ind_ix2 {
                Some(id) => {
                    assert!(id < inds.v().len() as u32);
                    id
                }
                None => {
                    let s = sample2.ind_name2.as_ref().unwrap();
                    inds.m()[s] as u32
                }
            };
            eprintln!(
                "Samples\n\tsample1:\t{}\t{}\n\tsample2:\t{}\t{}",
                inds.v()[id1 as usize],
                id1,
                inds.v()[id2 as usize],
                id2
            );

            if id1 < id2 {
                std::mem::swap(&mut id1, &mut id2);
            }

            let first = ibd1
                .as_slice()
                .partition_point(|x| x.individual_pair() < (id1, id2));
            let last = ibd1
                .as_slice()
                .partition_point(|x| x.individual_pair() <= (id1, id2));
            let v1 = &ibd1.as_slice()[first..last];

            let first = ibd2
                .as_slice()
                .partition_point(|x| x.individual_pair() < (id1, id2));
            let last = ibd2
                .as_slice()
                .partition_point(|x| x.individual_pair() <= (id1, id2));
            let v2 = &ibd2.as_slice()[first..last];

            let out = PathBuf::from(out).with_extension("svg");

            plot_svg(v1, v2, out, &ginfo).unwrap();
        }

        _ => {
            eprintln!("\nUse '-h  or [subcommand] -h' to show help message");
            std::process::exit(-1);
        }
    }
}

fn write_pair_total(
    pairtotal1: impl Iterator<Item = f32>,
    pairtotal2: impl Iterator<Item = f32>,
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

fn plot_svg(
    v1: &[IbdSeg],
    v2: &[IbdSeg],
    out: impl AsRef<Path>,
    ginfo: &GenomeInfo,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::{prelude::*, style::text_anchor::*};
    let root_area = SVGBackend::new(out.as_ref(), (1024, 768)).into_drawing_area();

    root_area.fill(&WHITE)?;

    let nchrom = ginfo.chromnames.len() as f32;
    let chrsz_max = *ginfo.chromsize.iter().max().unwrap() as f32;
    let mut cc = ChartBuilder::on(&root_area)
        .margin(40)
        .set_left_and_bottom_label_area_size(50)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .caption("IBD shared by a pair of samples", ("sans-serif", 30))
        .build_cartesian_2d(-1.0f32..chrsz_max, -1.0f32..nchrom)?;

    cc.configure_mesh()
        .x_labels(20)
        .y_labels(0)
        .disable_mesh()
        .x_desc("Position")
        .y_desc("Chromosome")
        .axis_desc_style(("sans-serif", 20))
        .x_label_formatter(&|v| format!("{:.0}", v))
        .y_label_formatter(&|v| format!("{:.0}", v))
        .draw()?;

    // draw IBD segments
    let mut points = vec![];
    for (iv, v) in [v1, v2].iter().enumerate() {
        eprintln!("Set{} size: {}", iv, v.len());
        for (iseg, seg) in v.iter().enumerate() {
            let (direction, color, label) = match iv {
                0 => (1.0f32, &BLUE, "Set1"),
                1 => (-1.0f32, &RED, "Set2"),
                _ => panic!(),
            };
            points.clear();
            let (_, m, _, n) = seg.haplotype_pair();
            let (chrid, chrname, s) = ginfo.to_chr_pos(seg.s);
            let e = seg.e - seg.s + s;
            eprintln!("\t{}\t{}\t{}", chrname, s + 1, e + 1);

            let y = chrid as f32 + (m + n + 1) as f32 * 0.03 * direction;
            points.push((s as f32, y));
            points.push((e as f32, y));

            let ls = LineSeries::new(
                points.iter().map(|x| (x.0, x.1)),
                color.clone().stroke_width(2),
            );

            let s = cc.draw_series(ls)?;
            let legstyle = |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color.clone());
            if iseg == 0 {
                s.label(label).legend(legstyle);
            }
        }
    }

    let right_center = Pos::new(HPos::Right, VPos::Center);

    // customize y tick labels
    cc.draw_series(PointSeries::of_element(
        (0.0f32..nchrom).step(1.0).values().map(|y| (0.0, y)),
        5,
        ShapeStyle::from(&BLACK).filled(),
        &|coord, size, style| {
            let ts = TextStyle {
                pos: right_center,
                font: ("sans-serif", 12).into(),
                color: RGBAColor(style.color.0, style.color.1, style.color.2, 0.8)
                    .to_backend_color(),
            };
            let chrid = coord.1 as usize;
            let chrname = &ginfo.chromnames[chrid];
            EmptyElement::at(coord)
                + PathElement::new(vec![(-size, 0), (0, 0)], style)
                + Text::new(format!("{}", chrname), (-size - 5, 0), ts)
        },
    ))?;

    cc.configure_series_labels().border_style(&BLACK).draw()?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root_area.present().expect("Unable to write result to file");
    println!(
        "Result has been saved to {}",
        out.as_ref().to_string_lossy()
    );
    Ok(())
}
