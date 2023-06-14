use clap::{Parser, Subcommand};
use ishare::{
    genome, gmap,
    indiv::*,
    share::ibd::{ibdset::*, overlap},
};
use std::path::PathBuf;

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
        // IBD directory 1
        ibd1_dir: PathBuf,
        // IBD directory 2
        ibd2_dir: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Some(Commands::Compare {
            genome_info,
            samples,
            out,
            ibd1_dir,
            ibd2_dir,
        }) => {
            // files
            let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
            let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
            let inds = Individuals::from_txt_file(&samples);
            let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds);
            ibd1.read_hapibd_dir(&ibd1_dir);
            let mut ibd2 = IbdSet::new(&gmap, &ginfo, &inds);
            ibd2.read_hapibd_dir(&ibd2_dir);

            ibd1.sort_by_samples();
            ibd2.sort_by_samples();

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
        _ => {
            println!("\nUse '-h  or [subcommand] -h' to show help message");
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
