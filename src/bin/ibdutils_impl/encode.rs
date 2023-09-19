use super::args::*;
use env_logger;
use ishare::{
    container::{intervals::Intervals, intervaltree::IntervalTree},
    genome::GenomeInfo,
    gmap::{self, GeneticMap},
    indiv::*,
    share::ibd::{ibdseg::IbdSeg, ibdset::*},
};
use log::*;
use std::path::PathBuf;

pub fn main_encode(args: &Commands) {
    if let Commands::Encode {
        genome_info,
        sample_lst,
        fmt,
        ibd_dir,
        outlier_upper,
        position_lst,
        min_snp,
        min_cm,
        out_prefix,
    } = args
    {
        env_logger::Builder::new()
            .filter(None, log::LevelFilter::Info)
            .format_module_path(false)
            .init();
        info!("read genome toml file");
        let ginfo = GenomeInfo::from_toml_file(&genome_info);
        info!("read genetic map files");
        let gmap = gmap::GeneticMap::from_genome_info(&ginfo);
        info!("read samples list file");
        let (inds, _inds_opt) = Individuals::from_txt_file(&sample_lst);

        // read ibd into memory
        let ibd = read_ibd(&ginfo, &gmap, &inds, ibd_dir, fmt);
        info!("before encoding no. ibd: {}", ibd.as_slice().len());

        // read positions
        let positions = read_positions(&ginfo, position_lst);

        // count snps per cM
        //   calculate low snp regions
        let (low_snp_regions, bounds, cnts) =
            get_low_snp_regions(&ginfo, &gmap, positions, min_snp);

        // get coverage
        let step_cm = 0.25;
        let (sp, cov) = get_coverage(&ibd, &gmap, step_cm);

        // calculate exetreme IBD sharing region
        let extreme_regions = get_exetreme_cov_regions(&sp, &cov, &ginfo, outlier_upper);

        // combined regions
        let mut region_to_remove = Intervals::new();
        region_to_remove.extend_from_iter(low_snp_regions.iter().map(|x| (x.start, x.end)));
        region_to_remove.extend_from_iter(extreme_regions.iter().map(|x| (x.start, x.end)));
        region_to_remove.merge();

        // complement
        let mut region_to_keep = region_to_remove.clone();
        region_to_keep.complement(0, ginfo.get_total_len_bp());

        // filtered IBD
        let mut v = cut_ibd(&ibd, &region_to_keep, *min_cm);

        // sort by haplotype
        v.sort();
        info!("after encoding no. ibd: {}", v.len());

        // write binary
        write_ibdseg_vec(&v, out_prefix);

        //    histogram: txt
        write_histogram(&sp, &cov, out_prefix, &ginfo);

        // write removed regions
        write_removed_region(&region_to_remove, out_prefix, &ginfo);
        // write snp density
        write_snp_counts(&bounds, &cnts, out_prefix, &ginfo);
    }
}
fn read_ibd<'a>(
    ginfo: &'a GenomeInfo,
    gmap: &'a GeneticMap,
    inds: &'a Individuals,
    ibd_dir: &PathBuf,
    fmt: &String,
) -> IbdSet<'a> {
    let mut ibd = IbdSet::new(&gmap, &ginfo, &inds);

    info!("read ibd list file");
    if fmt.as_str() == "hapibd" {
        ibd.read_hapibd_dir(ibd_dir);
    } else if fmt.as_str() == "tskibd" {
        ibd.read_tskibd_dir(ibd_dir);
    } else if fmt.as_str() == "hmmibd" {
        ibd.read_hmmibd_dir(ibd_dir);
    } else {
        panic!("format {} is not supported.", fmt);
    }
    ibd
}
fn read_positions(ginfo: &GenomeInfo, position_lst: &PathBuf) -> Vec<u32> {
    let s = std::fs::read_to_string(position_lst).expect("cannot read position list file");
    let mut positions = vec![];
    for lines in s.trim().split("\n") {
        let mut iter = lines.split("\t");
        let chrname = iter.next().expect("error getting chrname");
        let chr_pos: u32 = iter
            .next()
            .expect("error getting chr pos")
            .parse()
            .expect("can not parse chr_pos as u32");
        let chrid = ginfo
            .idx
            .get(chrname)
            .expect("chrname does not match genome toml file");
        let gw_pos = ginfo.to_gw_pos(*chrid, chr_pos);
        positions.push(gw_pos);
    }
    positions.sort();
    positions
}
fn get_low_snp_regions(
    ginfo: &GenomeInfo,
    gmap: &GeneticMap,
    positions: Vec<u32>,
    min_snp: &u32,
) -> (Intervals<u32>, Vec<u32>, Vec<usize>) {
    // boundaries in cM
    let boundaries = {
        let gwsize_cm = gmap.get_size_cm();
        let gwsize_bp = ginfo.get_total_len_bp();
        let mut x = 0.0f32;
        let mut v = vec![];
        while x < gwsize_cm {
            let bp = gmap.get_bp(x);
            v.push(bp);
            x += 1.0;
        }
        if *v.last().unwrap() < gwsize_bp {
            v.push(gwsize_bp);
        }
        v
    };

    // boundaries in bp
    let counts = {
        let mut counts = vec![0usize; boundaries.len() - 1];
        for p in positions.iter() {
            let idx = boundaries.partition_point(|x| x <= p) - 1;
            counts[idx] += 1;
        }
        counts
    };

    // low snp regions
    let low_snp_regions = {
        let mut regions = Intervals::new();
        for (i, count) in counts.iter().enumerate() {
            if *count < *min_snp as usize {
                let s = boundaries[i];
                let e = boundaries[i + 1];
                regions.push(s..e);
            }
        }
        regions
    };

    (low_snp_regions, boundaries, counts)
}

fn get_coverage(ibd: &IbdSet, gmap: &GeneticMap, step_cm: f32) -> (Vec<u32>, Vec<usize>) {
    // get sampling points
    let sp = {
        let mut sp = vec![];
        let mut x = step_cm / 2.0;
        let genome_size_cm = gmap.get_size_cm();
        // let genome_size_bp = ginfo.get_total_len_bp();
        while x < genome_size_cm {
            let gw_pos = gmap.get_bp(x);
            sp.push(gw_pos);
            x += step_cm;
        }
        sp
    };

    // coverage
    let cov = {
        let mut cov = vec![0; sp.len()];
        for seg in ibd.as_slice() {
            let i = sp.partition_point(|x| *x < seg.s) - 1;
            let j = sp.partition_point(|x| *x <= seg.e);
            cov[i..j].iter_mut().for_each(|x| *x += 1);
        }
        cov
    };
    (sp, cov)
}
fn get_exetreme_cov_regions(
    sp: &Vec<u32>,
    cov: &Vec<usize>,
    ginfo: &GenomeInfo,
    outlier_upper: &f64,
) -> Intervals<u32> {
    // expand data
    #[derive(Ord, PartialOrd, PartialEq, Eq)]
    struct Cov {
        chrid: usize,
        sp: u32,
        cv: usize,
        extreme: bool,
    }
    let mut data = Vec::with_capacity(sp.len());
    for (sp, cv) in sp.iter().zip(cov.iter()) {
        let (chrid, _, _) = ginfo.to_chr_pos(*sp);
        data.push(Cov {
            chrid,
            sp: *sp,
            cv: *cv,
            extreme: false,
        });
    }
    data.sort();

    use slice_group_by::GroupByMut;

    // mark extreme sp

    let mut buffer = vec![];
    for grp in data.as_mut_slice().linear_group_by_key_mut(|x| x.chrid) {
        buffer.clear();
        grp.iter().for_each(|x| buffer.push(x.cv));
        buffer.sort();
        // %3 trimmed
        let n = buffer.len() * 3 / 100;
        let buffer = &buffer[n..(buffer.len() - n)];
        // mean and std
        let n = buffer.len() as f64;
        let mean = buffer.iter().map(|x| *x as f64).sum::<f64>() / n;
        let mut std = buffer
            .iter()
            .map(|x| (*x as f64 - mean) * (*x as f64 - mean))
            .sum::<f64>()
            / n;
        std = std.sqrt();
        let threshold = mean + std * outlier_upper;

        grp.iter_mut().for_each(|x| {
            if x.cv as f64 > threshold {
                x.extreme = true;
            }
        })
    }

    // merge via interval (idx)
    let mut idx_intrvl = Intervals::new();
    for (i, d) in data.iter().enumerate() {
        if d.extreme {
            idx_intrvl.push(i..(i + 1));
        }
    }
    idx_intrvl.merge();

    // convert idx interval to positon interval
    let mut pos_intrvl = Intervals::new();
    for r in idx_intrvl.iter() {
        let s = r.start;
        let s_pos = data[s].sp;
        let e = r.end;
        let e_pos = if e < data.len() {
            data[e].sp
        } else {
            data.last().unwrap().sp
        };

        pos_intrvl.push(s_pos..e_pos);
    }

    pos_intrvl
}
// cut IBD by regions
fn cut_ibd(ibd: &IbdSet, region_to_keep: &Intervals<u32>, min_cm: f32) -> Vec<IbdSeg> {
    let mut out = Vec::<IbdSeg>::new();
    let tree = IntervalTree::from_iter(region_to_keep.iter().map(|x| ((x.start..x.end), ())));
    let gmap = ibd.get_gmap();

    for seg in ibd.iter() {
        for region in tree.query(seg.s..seg.e) {
            let mut s = region.range.start;
            let mut e = region.range.end;
            if s < seg.s {
                s = seg.s;
            }
            if e > seg.e {
                e = seg.e;
            }
            let mut new_seg = seg.clone();
            if gmap.get_cm_len(s, e) < min_cm {
                continue;
            }
            new_seg.s = s;
            new_seg.e = e;

            out.push(new_seg);
        }
    }
    out
}

fn write_ibdseg_vec(v: &Vec<IbdSeg>, out_prefix: &PathBuf) {
    use std::io::Write;
    let s = out_prefix.to_str().unwrap();
    let out = format!("{}.eibd", s);
    let mut file = std::fs::File::create(out.clone())
        .map(std::io::BufWriter::new)
        .expect(&format!("cannot create file {}", out));
    let sz = v.len() as u64;
    file.write(&sz.to_le_bytes()).unwrap();

    for seg in v.iter() {
        file.write(&seg.i.to_le_bytes()).unwrap();
        file.write(&seg.j.to_le_bytes()).unwrap();
        file.write(&seg.s.to_le_bytes()).unwrap();
        file.write(&seg.e.to_le_bytes()).unwrap();
    }
}

fn write_histogram(sp: &Vec<u32>, cov: &Vec<usize>, out_prefix: &PathBuf, ginfo: &GenomeInfo) {
    let s = out_prefix.to_str().unwrap();
    let out = format!("{}_hist.tsv", s);
    let mut file =
        std::fs::File::create(out.clone()).expect(&format!("cannot create file {}", out));
    use std::io::Write;
    for (s, c) in sp.iter().zip(cov.iter()) {
        let (_, name, pos) = ginfo.to_chr_pos(*s);
        write!(file, "{}\t{}\t{}\n", name, pos, c).unwrap();
    }
}
fn write_snp_counts(
    bounds: &Vec<u32>,
    counts: &Vec<usize>,
    out_prefix: &PathBuf,
    ginfo: &GenomeInfo,
) {
    let s = out_prefix.to_str().unwrap();
    let out = format!("{}_snp_counts.tsv", s);
    let mut file =
        std::fs::File::create(out.clone()).expect(&format!("cannot create file {}", out));
    use std::io::Write;
    for (s, c) in bounds.iter().zip(counts.iter()) {
        let (_, name, pos) = ginfo.to_chr_pos(*s);
        write!(file, "{}\t{}\t{}\n", name, pos, c).unwrap();
    }
}
fn write_removed_region(
    region_to_remove: &Intervals<u32>,
    out_prefix: &PathBuf,
    ginfo: &GenomeInfo,
) {
    let s = out_prefix.to_str().unwrap();
    let out = format!("{}_rm_regions.tsv", s);
    let mut file =
        std::fs::File::create(out.clone()).expect(&format!("cannot create file {}", out));
    use std::io::Write;
    for r in region_to_remove.iter() {
        let (_, name, s) = ginfo.to_chr_pos(r.start);
        let e = s + (r.end - r.start);
        write!(file, "{}\t{}\t{}\n", name, s, e).unwrap();
    }
}
