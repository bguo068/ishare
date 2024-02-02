use super::super::Commands;

use ishare::{
    genome, gmap,
    indiv::*,
    share::ibd::{coverage::CovCounter, ibdset::*},
};

pub fn main_coverage(args: &Commands) {
    if let Commands::Coverage {
        genome_info,
        sample_lst,
        fmt,
        ibd_dir,
        min_cm,
        start_cm,
        step_cm,
        prevent_flatten,
        out,
    } = args
    {
        let ginfo = genome::GenomeInfo::from_toml_file(&genome_info);
        let gmap = gmap::GeneticMap::from_genome_info(&ginfo);

        let (inds1, inds1_opt) = Individuals::from_txt_file(&sample_lst);
        let mut ibd1 = IbdSet::new(&gmap, &ginfo, &inds1);

        for (((ibd, fmt), dir), inds_opt) in [&mut ibd1]
            .iter_mut()
            .zip([fmt])
            .zip([ibd_dir])
            .zip([&inds1_opt])
        {
            if fmt.as_str() == "hapibd" {
                ibd.read_hapibd_dir(dir);
                if (*prevent_flatten) && (ibd.get_ploidy_status() == IbdSetPloidyStatus::Diploid) {
                    ibd.merge();
                }
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

            // remove short IBD segment to prevent unfair genome-wide overlapping analysis
            // for instance, original hmmIBD output all short IBD segments without filtering;
            // while other ibd caller such hap-IBD by default only output cm >= 2.0.
            // Directly filtering IBD segments by cm length take additonal bp to cM mapping.
            // Adding it here seems to be convenient
            ibd.filter_segments_by_min_cm(*min_cm);

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

        let total_size = gmap.get_size_cm();
        let mut sampling_points_cm: Vec<f32> = vec![];
        let mut cm_sp = *start_cm as f32;
        while cm_sp < total_size {
            sampling_points_cm.push(cm_sp);
            cm_sp += *step_cm as f32;
        }
        let sampling_points_bp: Vec<u32> = sampling_points_cm
            .iter()
            .map(|cm| gmap.get_bp(*cm as f32))
            .collect();

        let rngs = sampling_points_bp.iter().map(|bp| (*bp)..(*bp + 1));
        let mut counter = CovCounter::new(rngs);

        eprintln!("calculating coverage");
        for seg in ibd1.as_slice() {
            counter.count_over_interval(&(seg.s..seg.e));
        }

        let mut file = std::fs::File::create(out)
            .map(std::io::BufWriter::new)
            .unwrap();

        use std::io::Write;
        write!(file, "Chrom,Pos,Cm,GwPos,GwCm,Coverage\n").unwrap();
        for (gw_bp, _, coverage) in counter.iter_sorted_start_end_count() {
            let gw_cm = gmap.get_cm(gw_bp);
            let (idx, chr_str, chr_pos) = ginfo.to_chr_pos(gw_bp);
            let gwstart_bp = ginfo.gwstarts[idx];
            let gwstart_cm = gmap.get_cm(gwstart_bp);
            let chr_cm = gw_cm - gwstart_cm;

            write!(
                file,
                "{chr_str},{chr_pos},{chr_cm},{gw_bp},{gw_cm},{coverage}\n"
            )
            .unwrap();
        }
    }
}
