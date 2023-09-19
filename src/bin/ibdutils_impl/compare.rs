use super::utils::*;
use ishare::{
    genome, gmap,
    indiv::*,
    share::ibd::{ibdset::*, overlap},
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
            let ignore_hap = if *use_hap_overlap { false } else { true };
            // overlapping analysis
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
            let res = oa.analzyze(Some(&[3.0f32, 4.0, 6.0, 10.0, 18.0]));
            res.to_csv(out);
        }
        {
            // total IBD analysis
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

            write_pair_total(it1, it2, out.with_extension("totibd"));
        }
    }
}
