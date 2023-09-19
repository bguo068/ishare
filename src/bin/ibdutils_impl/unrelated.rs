use super::args::*;
use ahash::{HashMap, HashMapExt, HashSet, HashSetExt};
use env_logger;
use ishare::{
    genome::GenomeInfo,
    gmap::{self, GeneticMap},
    indiv::*,
    share::ibd::ibdset::*,
};
use log::*;
use std::path::PathBuf;

pub fn main_unrelated(args: &Commands) {
    if let Commands::GetUnrelated {
        genome_info,
        sample_lst,
        fmt,
        ibd_dir,
        theshold,
        out,
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
        // read ibd
        let ibd = read_ibd(&ginfo, &gmap, &inds, ibd_dir, fmt);
        info!("total IBD segment counts: {}", ibd.as_slice().len());

        let related_paris = get_related_pairs(&ibd, *theshold);
        let samples_to_rm = get_samples_to_remove(related_paris);
        info!(
            "remove {} samples out of {}",
            samples_to_rm.len(),
            inds.v().len()
        );
        write_sample_to_keep(out, &inds, samples_to_rm);
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
    ibd.infer_ploidy();
    // this will also sort by samples
    info!("sort and merge ibd");
    ibd.merge();
    ibd
}

fn get_related_pairs(ibd: &IbdSet, threshold: f64) -> Vec<(u32, u32)> {
    let gmap = ibd.get_gmap();
    let gsize = gmap.get_size_cm();
    let min_totibd_related = gsize * threshold as f32;
    info!(
        "related pair total ibd theshold in cM: {:.3}",
        min_totibd_related
    );

    let mut v = vec![];
    for blk in IbdSetBlockIter::new(ibd, true) {
        let totibd: f32 = blk.iter().map(|seg| seg.get_seg_len_cm(&gmap)).sum();
        if totibd <= min_totibd_related {
            continue;
        }
        let related_pair = blk.first().unwrap().individual_pair();
        v.push(related_pair);
    }
    info!("found related pairs:: {}", v.len());
    v
}

fn get_samples_to_remove(mut related_pairs: Vec<(u32, u32)>) -> HashSet<u32> {
    let mut hashset = HashSet::new();
    let mut degree = HashMap::<u32, u32>::with_capacity(related_pairs.len() * 2);

    while related_pairs.len() > 0 {
        calc_degrees(&related_pairs, &mut degree);
        let sample_with_max_degree = get_sample_with_largest_degree(&degree);
        remove_pairs_with_sample(&mut related_pairs, sample_with_max_degree);
        hashset.insert(sample_with_max_degree);
    }

    hashset
}

fn calc_degrees(related_paris: &Vec<(u32, u32)>, degree: &mut HashMap<u32, u32>) {
    // clear degree
    degree.clear();
    for (i, j) in related_paris.iter() {
        let cnt = degree.entry(*i).or_insert(0);
        *cnt += 1;
        let cnt = degree.entry(*j).or_insert(0);
        *cnt += 1;
    }
}

fn get_sample_with_largest_degree(degree: &HashMap<u32, u32>) -> u32 {
    let mut sample = 0;
    let mut max_deg = 0;
    for (k, v) in degree {
        if *v > max_deg {
            max_deg = *v;
            sample = *k;
        }
    }
    sample
}

fn remove_pairs_with_sample(related_paris: &mut Vec<(u32, u32)>, target_sample: u32) {
    for (i, j) in related_paris.iter_mut() {
        if (*i == target_sample) || (*j == target_sample) {
            // mark for removal;
            *i = u32::MAX;
            *j = u32::MAX;
        }
    }
    // consolidate
    related_paris.retain(|x| !((x.0 == u32::MAX) & (x.1 == u32::MAX)));
}

fn write_sample_to_keep(out: &PathBuf, inds: &Individuals, samples_to_rm: HashSet<u32>) {
    use std::io::Write;
    let mut file = std::fs::File::create(out)
        .map(std::io::BufWriter::new)
        .unwrap();
    for (i, name) in inds.v().iter().enumerate() {
        if !samples_to_rm.contains(&(i as u32)) {
            write!(file, "{}\n", name).expect("faile to write sample name to file");
        }
    }
}
