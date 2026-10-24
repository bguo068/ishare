use super::{GtencodeError, Result};
use error_stack::*;
use skato::*;

use super::Commands;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, site::Sites};

pub fn main_skato(args: &Commands) -> Result<()> {
    if let Commands::Skato {
        rec,
        genomes: _,
        window,
        step,
    } = args
    {
        // read files
        let mut records =
            GenotypeRecords::from_parquet_file(rec).change_context(GtencodeError::Input)?;

        let ind_file = rec.with_extension("ind");
        let inds =
            Individuals::from_parquet_file(&ind_file).change_context(GtencodeError::Input)?;

        let sites_file = rec.with_extension("sit");
        let sites = Sites::from_parquet_file(&sites_file).change_context(GtencodeError::Input)?;

        println!("number of individuals: {}", inds.v().len());
        println!("number of sites      : {}", sites.len());

        // check sort status
        println!("check sort status of genotype records");
        if !records.is_sorted_by_position_genome_allele_slow() {
            records.assume_unsorted();
            println!("sorting gennotype");
            records
                .sort_by_position_genome_allele()
                .change_context(GtencodeError::Library)?;
        }

        // data slices
        let all_records = records.records();
        let sitepos = sites.get_gw_pos_slice();
        assert!(sitepos.is_sorted());

        let it1 = (0..sitepos.len()).step_by(*step);
        let it2 = (0..sitepos.len()).skip(*window).step_by(*step);

        let mut total_steps = (sitepos.len() - *window) / (*step);
        if ((sitepos.len() - *window) % (*step)) != 0 {
            total_steps += 1;
        }
        // others --------------------
        let nsams = inds.v().len();
        println!("\trun logistic regression");
        // make a fake x matrix
        let mut x = Array2::<f64>::from_elem([nsams, 1], 1.0);
        x.iter_mut()
            .enumerate()
            .for_each(|(i, x)| *x = (i / 2) as f64);

        // make a fake y matrix
        let y = Array1::<f64>::from_iter((0..nsams).map(|n| (n % 2) as f64));
        // MU
        let mu = run_logistic_regression_b(&x, &y);
        // PI
        // println!("get PI");
        let pi: Array1<f64> = mu.clone().map(|mu_i| mu_i * (1.0 - mu_i));
        let rho_vec = vec![0.1, 0.3, 0.5, 0.7, 0.9];

        // iterator windows
        for (istep, (min_idx, max_idx)) in it1.zip(it2).enumerate() {
            let min_pos = sitepos[min_idx];
            let max_pos = sitepos[max_idx];
            let span_pos = max_pos - min_pos;
            println!(
                "istep = {istep}, percentage = {:.4} ======================",
                istep as f64 / total_steps as f64
            );
            println!("\tsite min_pos = {min_pos}, site max_pos = {max_pos}, span = {span_pos}");
            println!("\tnumber of selected sites  : {}", *window);

            let sel_records = {
                let s = all_records.partition_point(|x| x.get_position() < min_pos);
                let e = all_records.partition_point(|x| x.get_position() < max_pos);
                // dbg!(s, e, min_pos, max_pos);
                &all_records[s..e]
            };
            if sel_records.is_empty() {
                continue;
            }
            println!("\tnumber of selected records: {}", sel_records.len());
            let nsites = *window;

            {
                println!("\tprepare matrix for call skat-o test");
                // build matrix
                let mut z = Array2::<f64>::zeros([nsams, nsites]);
                // get genotype data Z
                println!("get genotype data Z");
                let mut target_pos_idx = 0;
                for r in sel_records {
                    let gw_pos = r.get_position();
                    let (global_idx, _global_idx_plus) = sites.get_idx_by_position(gw_pos);
                    // dbg!(global_idx, global_idx_plus);
                    // idx within selected sites
                    let pos_idx = global_idx - min_idx;
                    let sam_idx = r.get_genome() / 2;
                    // assert!(sam_idx < nsams as u32);
                    // dbg!(pos_idx, nsites, sitepos.len());
                    assert!(pos_idx < nsites,);
                    z[[sam_idx as usize, pos_idx]] += 1.0;
                    if (pos_idx != 0) && (target_pos_idx == 0) {
                        target_pos_idx = pos_idx;
                    }
                }
                // MAF
                println!("get MAF");
                let maf = get_maf(&z);
                // WEIGHTS
                println!("get WEIGHTS");
                let mut weights = Vec::new();
                get_beta_weights(&maf, 1.0, 25.0, &mut weights);
                let weights: Array1<f64> = weights.into();
                // X1
                println!("get X1");
                let x1 = make_mat_x1(&x);

                println!("\trun Skat-O test");
                let mut calculator = SkatOCalculator::new(&x1, z, &weights, &y, &mu, &pi, &rho_vec);
                println!("\tmultiply_weight_to_z");
                calculator.multiply_weight_to_z();
                println!("\tmultiply_p_to_z");
                calculator.multiply_p_to_z();
                // do this step before z multiply l (which is called in get_pmin_and_quantiles)
                println!("\tget_t_params");
                let t_params = calculator.get_t_params();
                println!("\tget_pmin_and_quantiles");
                let (_t, quantiles) = calculator.get_pmin_and_quantiles();
                println!("\tget_t_p_value");
                let t_pval = calculator.get_t_p_value(&quantiles, &t_params);

                println!("t_pval = {t_pval}");
            }
        }
    }
    Ok(())
}
