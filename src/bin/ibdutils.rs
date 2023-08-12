use clap::{Args, Parser, Subcommand};
use ishare::{
    genome::{self, GenomeInfo},
    gmap,
    indiv::*,
    share::ibd::{ibdseg::IbdSeg, ibdset::*, overlap},
};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name="ibdutils", author, version, about, long_about=None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Encode genotype data from VCF file to tables or matrix (-m)
    Compare {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        /// Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst1: PathBuf,
        /// Path to sample list file for ibd set2
        #[arg(short = 'S', long, required = true)]
        sample_lst2: PathBuf,
        /// fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt1: String,
        /// fmt of ibd set2, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'F', long, required = true)]
        fmt2: String,
        /// IBD directory 1
        #[arg(short = 'i', long, required = true)]
        ibd1_dir: PathBuf,
        /// IBD directory 2
        #[arg(short = 'I', long, required = true)]
        ibd2_dir: PathBuf,
        /// Use (Do not Ignore) haplotype information (i.e. not flattening before overlapping analysis) if true
        #[arg(short = 'N', long, default_value_t = false)]
        use_hap_info: bool,
        /// Path to sample list file
        #[arg(short = 'o', long, default_value = "ibd_cmp_res.txt")]
        out: PathBuf,
    },

    PlotIBD {
        /// Path to genome info toml file (input)
        #[arg(short = 'g', long, default_value = "genome.toml")]
        genome_info: PathBuf,
        // Path to sample list file for ibd set1
        #[arg(short = 's', long, required = true)]
        sample_lst1: PathBuf,
        // Path to sample list file for ibd set2
        #[arg(short = 'S', long, required = true)]
        sample_lst2: PathBuf,
        // fmt of ibd set1, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'f', long, required = true)]
        fmt1: String,
        // fmt of ibd set2, supported format 'hapibd', 'tskibd' and 'hmmibd'
        #[arg(short = 'F', long, default_value = "hapibd")]
        fmt2: String,
        // IBD directory 1
        #[arg(short = 'i', long, required = true)]
        ibd1_dir: PathBuf,
        // IBD directory 2
        #[arg(short = 'I', long, required = true)]
        ibd2_dir: PathBuf,
        // Path to sample list file
        #[arg(short = 'o', long, default_value = "ibd_cmp_res.txt")]
        out: PathBuf,

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
            sample_lst1,
            sample_lst2,
            fmt1,
            fmt2,
            ibd1_dir,
            ibd2_dir,
            use_hap_info,
            out,
        }) => {
            // files
            let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
            let gmap = gmap::GeneticMap::from_genome_info(&ginfo);

            let (inds1, inds1_opt) = Individuals::from_txt_file(&sample_lst1);
            let (inds2, inds2_opt) = Individuals::from_txt_file(&sample_lst2);
            let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds1);
            let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds2);

            for (((ibd, fmt), dir), inds_opt) in [&mut ibd1, &mut ibd2]
                .iter_mut()
                .zip([fmt1, fmt2])
                .zip([ibd1_dir, ibd2_dir])
                .zip([&inds1_opt, &inds2_opt])
            {
                if fmt.as_str() == "hapibd" {
                    ibd.read_hapibd_dir(dir);
                    ibd.sort_by_samples();
                    ibd.infer_ploidy();
                } else if fmt.as_str() == "tskibd" {
                    ibd.read_tskibd_dir(dir);
                    ibd.sort_by_haplotypes();
                    ibd.infer_ploidy();
                } else if fmt.as_str() == "hmmibd" {
                    ibd.read_hmmibd_dir(dir);
                    ibd.sort_by_haplotypes();
                    ibd.infer_ploidy();
                } else {
                    panic!("format {} is not supported.", fmt);
                }
                match inds_opt.as_ref() {
                    Some((converter, ind_, PloidConvertDirection::Diploid2Haploid)) => {
                        ibd.covert_to_haploid(ind_, converter);
                    }
                    Some((converter, ind_, PloidConvertDirection::Haploid2Diploid)) => {
                        ibd.covert_to_het_diploid(ind_, converter);
                    }
                    None => {}
                }
            }

            let ignore_hap = if *use_hap_info { false } else { true };
            assert_eq!(ibd1.get_inds().v(), ibd2.get_inds().v());

            {
                // overlapping analysis
                let oa = overlap::IbdOverlapAnalyzer::new(&mut ibd1, &mut ibd2, ignore_hap);
                let res = oa.analzyze(Some(&[3.0f32, 4.0, 6.0, 10.0, 18.0]));
                res.to_csv(out);
            }
            {
                let inds = ibd1.get_inds();
                // total IBD analysis
                let mat1 = ibd1.get_gw_total_ibd_matrix(ignore_hap);
                let mat2 = ibd2.get_gw_total_ibd_matrix(ignore_hap);
                let mut n = inds.v().len();
                if !ignore_hap {
                    n *= 2; // n is the number of haplotypes if not ignore_hap
                }
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
            sample_lst1,
            sample_lst2,
            fmt1,
            fmt2,
            ibd1_dir,
            ibd2_dir,
            out,
            sample1,
            sample2,
        }) => {
            let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
            let gmap = gmap::GeneticMap::from_genome_info(&ginfo);

            let (inds1, inds1_opt) = Individuals::from_txt_file(&sample_lst1);
            let (inds2, inds2_opt) = Individuals::from_txt_file(&sample_lst2);
            let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds1);
            let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds2);

            for (((ibd, fmt), dir), inds_opt) in [&mut ibd1, &mut ibd2]
                .iter_mut()
                .zip([fmt1, fmt2])
                .zip([ibd1_dir, ibd2_dir])
                .zip([&inds1_opt, &inds2_opt])
            {
                if fmt.as_str() == "hapibd" {
                    ibd.read_hapibd_dir(dir);
                    ibd.sort_by_samples();
                    ibd.infer_ploidy();
                    println!("{:?}", ibd.as_slice().len());
                } else if fmt.as_str() == "tskibd" {
                    ibd.read_tskibd_dir(dir);
                    ibd.sort_by_haplotypes();
                    ibd.infer_ploidy();
                } else if fmt.as_str() == "hmmibd" {
                    ibd.read_hmmibd_dir(dir);
                    ibd.sort_by_haplotypes();
                    ibd.infer_ploidy();
                } else {
                    panic!("format {} is not supported.", fmt);
                }
                match inds_opt.as_ref() {
                    Some((converter, ind_, PloidConvertDirection::Diploid2Haploid)) => {
                        ibd.covert_to_haploid(ind_, converter);
                    }
                    Some((converter, ind_, PloidConvertDirection::Haploid2Diploid)) => {
                        ibd.covert_to_het_diploid(ind_, converter);
                    }
                    None => {}
                }
            }

            assert_eq!(ibd1.get_inds().v(), ibd2.get_inds().v());
            let inds = ibd1.get_inds();

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
