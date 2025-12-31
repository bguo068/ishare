use super::args::*;
use env_logger;
use ishare::{
    container::{intervals::Intervals, intervaltree::IntervalTree},
    genome::GenomeInfo,
    gmap::{self, GeneticMap},
    indiv::*,
    share::ibd::{ibdseg::*, ibdset::*},
};
use itertools::Itertools;
use log::*;
use std::sync::Arc;
use std::{
    backtrace::Backtrace,
    path::{Path, PathBuf},
};

use snafu::prelude::*;
#[derive(Debug, Snafu)]
pub enum Error {
    // #[snafu(transparent)]
    Indiv {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::indiv::Error,
    },
    // #[snafu(transparent)]
    Genome {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::genome::Error,
    },
    // #[snafu(transparent)]
    Gmap {
        // non leaf
        #[snafu(source(from(ishare::gmap::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<ishare::gmap::Error>,
    },
    // #[snafu(transparent)]
    Ibd {
        // non leaf
        #[snafu(source(from(ishare::share::ibd::Error, Box::new)))]
        #[snafu(backtrace)]
        source: Box<ishare::share::ibd::Error>,
    },
    // #[snafu(transparent)]
    Container {
        // non leaf
        #[snafu(backtrace)]
        source: ishare::container::Error,
    },
    // local
    // #[snafu(transparent)]
    StdIo {
        // leaf
        source: std::io::Error,
        backtrace: Box<Option<Backtrace>>,
    },
    ReadPositions {
        // leaf
        msg: Box<String>,
        backtrace: Box<Option<Backtrace>>,
    },
    LowSnpWindowEmpty {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
    ZeroSamplingPoint {
        // leaf
        backtrace: Box<Option<Backtrace>>,
    },
}
type Result<T> = std::result::Result<T, Error>;

pub fn main_encode(args: &Commands) -> Result<()> {
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
        let ginfo = Arc::new(GenomeInfo::from_toml_file(genome_info).context(GenomeSnafu)?);
        info!("read genetic map files");
        let gmap = Arc::new(gmap::GeneticMap::from_genome_info(&ginfo).context(GmapSnafu)?);
        info!("read samples list file");
        let (inds, _inds_opt) = Individuals::from_txt_file(sample_lst).context(IndivSnafu)?;
        let inds = Arc::new(inds);

        // read ibd into memory
        let ibd = read_ibd(ginfo.clone(), gmap.clone(), inds, ibd_dir, fmt)?;
        info!("before encoding no. ibd: {}", ibd.as_slice().len());

        // read positions
        let positions = match position_lst.as_ref() {
            Some(position_lst) => read_positions(&ginfo, position_lst)?,
            None => get_uniq_positions_from_ibd(&ibd),
        };

        // count snps per cM
        //   calculate low snp regions
        let (low_snp_regions, bounds, cnts) =
            get_low_snp_regions(&ginfo, &gmap, positions, min_snp)?;

        // get coverage
        let step_cm = 0.25;
        let (sp, cov) = get_coverage(&ibd, &gmap, step_cm)?;

        // calculate exetreme IBD sharing region
        let extreme_regions = get_exetreme_cov_regions(&sp, &cov, &ginfo, outlier_upper)?;

        // combined regions
        let mut region_to_remove = Intervals::new();
        region_to_remove.extend_from_iter(low_snp_regions.iter().map(|x| (x.start, x.end)));
        region_to_remove.extend_from_iter(extreme_regions.iter().map(|x| (x.start, x.end)));
        region_to_remove.merge();

        // complement
        let mut region_to_keep = region_to_remove.clone();
        region_to_keep
            .complement(0, ginfo.get_total_len_bp())
            .context(ContainerSnafu)?;

        // filtered IBD
        let mut v = cut_ibd(&ibd, &region_to_keep, *min_cm);

        // sort by haplotype
        v.sort();
        info!("after encoding no. ibd: {}", v.len());

        // write binary
        let out = format!("{}.eibd", out_prefix.to_string_lossy());
        write_ibdseg_vec(&v, &out).context(IbdSnafu)?;

        //    histogram: txt
        write_histogram(&sp, &cov, out_prefix, &ginfo)?;

        // write removed regions
        write_removed_region(&region_to_remove, out_prefix, &ginfo)?;
        // write snp density
        write_snp_counts(&bounds, &cnts, out_prefix, &ginfo)?;
    }
    Ok(())
}
fn read_ibd(
    ginfo: Arc<GenomeInfo>,
    gmap: Arc<GeneticMap>,
    inds: Arc<Individuals>,
    ibd_dir: &PathBuf,
    fmt: &String,
) -> Result<IbdSet> {
    let mut ibd = IbdSet::new(gmap, ginfo, inds);

    info!("read ibd list file");
    if fmt.as_str() == "hapibd" {
        ibd.read_hapibd_dir(ibd_dir).context(IbdSnafu)?;
    } else if fmt.as_str() == "tskibd" {
        ibd.read_tskibd_dir(ibd_dir).context(IbdSnafu)?;
    } else if fmt.as_str() == "hmmibd" {
        ibd.read_hmmibd_dir(ibd_dir).context(IbdSnafu)?;
    } else {
        panic!("format {fmt} is not supported.");
    }
    Ok(ibd)
}
fn read_positions(ginfo: &GenomeInfo, position_lst: &PathBuf) -> Result<Vec<u32>> {
    let s = std::fs::read_to_string(position_lst).context(StdIoSnafu)?;
    let mut positions = vec![];
    for lines in s.trim().split("\n") {
        let mut iter = lines.split("\t");
        let chrname = iter.next().context(ReadPositionsSnafu {
            msg: "error getting chrname".to_owned(),
        })?;
        let chr_pos: u32 = iter
            .next()
            .context(ReadPositionsSnafu {
                msg: "error getting chr pos".to_owned(),
            })?
            .parse::<u32>()
            .ok()
            .context(ReadPositionsSnafu {
                msg: "can not parse chr_pos as u32".to_owned(),
            })?;
        let chrid = ginfo.idx.get(chrname).context(ReadPositionsSnafu {
            msg: "chrname does not match genome toml file".to_owned(),
        })?;
        let gw_pos = ginfo.to_gw_pos(*chrid, chr_pos);
        positions.push(gw_pos);
    }
    positions.sort();
    Ok(positions)
}
fn get_low_snp_regions(
    ginfo: &GenomeInfo,
    gmap: &GeneticMap,
    positions: Vec<u32>,
    min_snp: &u32,
) -> Result<(Intervals<u32>, Vec<u32>, Vec<usize>)> {
    // boundaries in cM
    let boundaries = {
        let gwsize_cm = gmap.get_size_cm().context(GmapSnafu)?;
        let gwsize_bp = ginfo.get_total_len_bp();
        let mut x = 0.0f32;
        let mut v = vec![];
        while x < gwsize_cm {
            let bp = gmap.get_bp(x);
            v.push(bp);
            x += 1.0;
        }
        if *v.last().context(LowSnpWindowEmptySnafu)? < gwsize_bp {
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

    Ok((low_snp_regions, boundaries, counts))
}

fn get_coverage(ibd: &IbdSet, gmap: &GeneticMap, step_cm: f32) -> Result<(Vec<u32>, Vec<usize>)> {
    // get sampling points
    let sp = {
        let mut sp = vec![];
        let mut x = step_cm / 2.0;
        let genome_size_cm = gmap.get_size_cm().context(GmapSnafu)?;
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
            let i = sp.partition_point(|x| *x < seg.s);
            let j = sp.partition_point(|x| *x <= seg.e);
            cov[i..j].iter_mut().for_each(|x| *x += 1);
        }
        cov
    };
    Ok((sp, cov))
}
fn get_exetreme_cov_regions(
    sp: &[u32],
    cov: &[usize],
    ginfo: &GenomeInfo,
    outlier_upper: &f64,
) -> Result<Intervals<u32>> {
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
            data.last().context(ZeroSamplingPointSnafu)?.sp
        };

        pos_intrvl.push(s_pos..e_pos);
    }

    Ok(pos_intrvl)
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
            let mut new_seg = *seg;
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

fn write_histogram(sp: &[u32], cov: &[usize], out_prefix: &Path, ginfo: &GenomeInfo) -> Result<()> {
    let s = out_prefix.to_string_lossy();
    let out = format!("{s}_hist.tsv");
    let mut file =
        std::fs::File::create(out.clone()).unwrap_or_else(|_| panic!("cannot create file {out}"));
    use std::io::Write;
    for (s, c) in sp.iter().zip(cov.iter()) {
        let (_, name, pos) = ginfo.to_chr_pos(*s);
        writeln!(file, "{name}\t{pos}\t{c}").context(StdIoSnafu)?;
    }
    Ok(())
}
fn write_snp_counts(
    bounds: &[u32],
    counts: &[usize],
    out_prefix: &Path,
    ginfo: &GenomeInfo,
) -> Result<()> {
    let s = out_prefix.to_string_lossy();
    let out = format!("{s}_snp_counts.tsv");
    let mut file =
        std::fs::File::create(out.clone()).unwrap_or_else(|_| panic!("cannot create file {out}"));
    use std::io::Write;
    for (s, c) in bounds.iter().zip(counts.iter()) {
        let (_, name, pos) = ginfo.to_chr_pos(*s);
        writeln!(file, "{name}\t{pos}\t{c}").context(StdIoSnafu)?;
    }
    Ok(())
}
fn write_removed_region(
    region_to_remove: &Intervals<u32>,
    out_prefix: &Path,
    ginfo: &GenomeInfo,
) -> Result<()> {
    let s = out_prefix.to_string_lossy();
    let out = format!("{s}_rm_regions.tsv");
    let mut file =
        std::fs::File::create(out.clone()).unwrap_or_else(|_| panic!("cannot create file {out}"));
    use std::io::Write;
    for r in region_to_remove.iter() {
        let (chrid_s, name_s, s) = ginfo.to_chr_pos(r.start);
        let (chrid_e, name_e, e) = ginfo.to_chr_pos(r.end);
        if chrid_s == chrid_e {
            writeln!(file, "{name_s}\t{s}\t{e}").context(StdIoSnafu)?;
        } else {
            let mut start;
            let mut end;
            let mut name;
            for chrid in chrid_s..=chrid_e {
                if chrid == chrid_s {
                    // first
                    start = s;
                    end = ginfo.chromsize[chrid_s] - 1;
                    name = name_s;
                } else if chrid == chrid_e {
                    // last
                    start = 0;
                    end = e;
                    name = name_e;
                } else {
                    // every full chromosomes in the middle
                    start = 0;
                    end = ginfo.chromsize[chrid] - 1;
                    name = ginfo.chromnames[chrid].as_str();
                }
                writeln!(file, "{name}\t{start}\t{end}").context(StdIoSnafu)?;
            }
        }
    }
    Ok(())
}
fn get_uniq_positions_from_ibd(ibd: &IbdSet) -> Vec<u32> {
    use ahash::{HashSet, HashSetExt};
    let mut set = HashSet::new();
    for seg in ibd.iter() {
        set.insert(seg.s);
        set.insert(seg.e);
    }
    let mut v = set.iter().copied().collect_vec();
    v.sort();
    v
}
