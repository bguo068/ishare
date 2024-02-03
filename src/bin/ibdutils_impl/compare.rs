use super::utils::*;
use ishare::{
    genome, gmap,
    indiv::*,
    share::ibd::{
        ibdset::*,
        overlap::{self, write_per_winddow_overlap_res, IbdOverlapResult},
    },
};

use super::super::Commands;
pub fn main_compare(args: &Commands) {
    if let Commands::Compare {
        genome_info,
        sample_lst1,
        sample_lst2,
        fmt1,
        fmt2,
        ibd1_dir,
        ibd2_dir,
        min_cm,
        length_bin_starts,
        window_size_bp,
        use_hap_overlap,
        use_hap_totibd,
        write_details,
        out,
    } = args
    {
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
                match *use_hap_overlap {
                    true => ibd.sort_by_haplotypes(),
                    false => ibd.sort_by_samples(),
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

        assert_eq!(ibd1.get_inds().v(), ibd2.get_inds().v());

        {
            // convert from string to vector of float
            let mut length_bin_starts: Vec<f32> = length_bin_starts
                .split(",")
                .map(|s| s.parse().unwrap())
                .collect();
            // sort
            length_bin_starts.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let ignore_hap = if *use_hap_overlap { false } else { true };

            // 1. overlapping analysis
            let prefix_for_details = match *write_details {
                true => Some(out),
                false => None,
            };
            let oa = overlap::IbdOverlapAnalyzer::new(
                &mut ibd1,
                &mut ibd2,
                ignore_hap,
                prefix_for_details,
            );

            // 1.1 overlapping analysis all together
            let res = oa.analzyze(Some(length_bin_starts.as_slice()), None);
            res.to_csv(out);

            // 1.2 overlapping analysis per window of each chromosome
            if let Some(window_size_bp) = window_size_bp {
                let nchrom = ginfo.gwstarts.len();
                let mut windows = Vec::<(u32, u32)>::new();
                let mut ov_res_windows = Vec::<IbdOverlapResult>::new();
                for idx in 0..nchrom {
                    let chr_gw_start = ginfo.gwstarts[idx];
                    let chr_gw_end = if idx == nchrom - 1 {
                        ginfo.chromsize[nchrom - 1] + ginfo.gwstarts[nchrom - 1]
                    } else {
                        ginfo.gwstarts[idx + 1]
                    };
                    let mut winstart = chr_gw_start;
                    while winstart < chr_gw_end {
                        let mut winend = winstart + window_size_bp - 1;
                        if winend > chr_gw_end {
                            winend = chr_gw_end - 1;
                        }
                        windows.push((winstart, winend));
                        let ov_res = oa
                            .analzyze(Some(length_bin_starts.as_slice()), Some((winstart, winend)));
                        ov_res_windows.push(ov_res);
                        winstart += window_size_bp;
                    }
                }
                // write to file
                let out_win_ov = out.with_extension("winovcsv");
                write_per_winddow_overlap_res(
                    &windows,
                    &ov_res_windows,
                    &ginfo,
                    &gmap,
                    &out_win_ov,
                );
            }
        }
        {
            // 2. total IBD analysis
            let ignore_hap = if *use_hap_totibd { false } else { true };

            // re-sort ibd
            match ignore_hap {
                true => {
                    for ibd in [&mut ibd1, &mut ibd2] {
                        if !ibd.is_sorted_by_samples() {
                            ibd.sort_by_samples();
                        }
                        ibd.infer_ploidy()
                    }
                }
                false => {
                    for ibd in [&mut ibd1, &mut ibd2] {
                        if !ibd.is_sorted_by_haplotypes() {
                            ibd.sort_by_haplotypes();
                        }
                        ibd.infer_ploidy()
                    }
                }
            }

            let inds = ibd1.get_inds();

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

            write_pair_total(
                it1,
                it2,
                "PairTotIbdA",
                "PairTotIbdB",
                out.with_extension("pairtotibdpq"),
            );
        }

        {
            // 3. population-level total ibd
            let binwidth = 0.05f32;
            let mut bincenters = vec![];
            let mut totals1 = vec![];
            let mut totals2 = vec![];

            for seg in ibd1.as_slice() {
                let cm = seg.get_seg_len_cm(&gmap);
                let idx = (cm / binwidth) as usize;
                if idx >= totals1.len() {
                    totals1.resize(idx + 1, 0.0f32);
                }
                totals1[idx] += cm;
            }

            for seg in ibd2.as_slice() {
                let cm = seg.get_seg_len_cm(&gmap);
                let idx = (cm / binwidth) as usize;
                if idx >= totals2.len() {
                    totals2.resize(idx + 1, 0.0f32);
                }
                totals2[idx] += cm;
            }

            // keep two vector same size
            let mut sz = totals1.len();
            if sz < totals2.len() {
                sz = totals2.len();
            }

            totals1.resize(sz, 0.0);
            totals2.resize(sz, 0.0);

            // make bin center vector
            for idx in 1..sz {
                let c = (idx as f32 + 0.5) * binwidth;
                bincenters.push(c);
            }

            // write csv file
            let mut f = std::fs::File::create(out.with_extension("poptotibdcsv")).unwrap();
            use std::io::Write;
            write!(f, "BinCenter,PopTotIbdA,PopTotIbdB\n").unwrap();
            for ((bc, tot1), tot2) in bincenters.iter().zip(totals1.iter()).zip(totals2.iter()) {
                write!(f, "{bc},{tot1},{tot2}\n").unwrap();
            }
        }
    }
}
