use super::super::Commands;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, site::Sites};

pub fn main_skato(args: &Commands) {
    if let Commands::Skato {
        rec,
        genomes: _,
        window,
        step,
    } = args
    {
        let mut records = GenotypeRecords::from_parquet_file(rec);
        let ind_file = rec.with_extension("ind");
        let inds = Individuals::from_parquet_file(&ind_file);
        let sites_file = rec.with_extension("sit");
        let sites = Sites::from_parquet_file(&sites_file);
        println!("number of individuals: {}", inds.v().len());
        println!("number of sites      : {}", sites.len());

        println!("sort by genotype records by position");
        records.sort_by_position();
        let all_records = records.records();
        let sitepos = sites.get_gw_pos_slice();

        let it1 = (0..sitepos.len()).step_by(*step as usize);
        let it2 = (0..sitepos.len())
            .skip(*window as usize)
            .step_by(*step as usize);

        let mut total_steps = (sitepos.len() - *window) / (*step);
        if ((sitepos.len() - *window) % (*step)) != 0 {
            total_steps += 1;
        }

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
                &all_records[s..e]
            };
            println!("\tnumber of selected records: {}", sel_records.len());
            let nsites = *window as usize;
            let nsams = inds.v().len();
            {
                println!("\tprepare matrix for call skat-o test");
                use skato_rs::*;
                // build matrix
                let mut z = Array2::<f64>::zeros([nsams, nsites]);
                // get genotype data Z
                let mut target_pos_idx = 0;
                for r in sel_records {
                    let gw_pos = r.get_position();
                    let global_idx = sites.get_idx_by_position(gw_pos).0;
                    // idx within selected sites
                    let pos_idx = global_idx - min_idx;
                    let sam_idx = r.get_genome() / 2;
                    z[[sam_idx as usize, pos_idx]] += 1.0;
                    if (pos_idx != 0) && (target_pos_idx == 0) {
                        target_pos_idx = pos_idx;
                    }
                }
                // make a fake x matrix
                let mut x = Array2::<f64>::from_elem([nsams, 1], 1.0);
                x.iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = (i / 2) as f64);

                // make a fake y matrix
                let y = z.column(target_pos_idx).to_owned();

                // others --------------------
                println!("\trun logistic regression");
                // MAF
                let maf = get_maf(&z);
                // WEIGHTS
                let mut weights = Vec::new();
                get_beta_weights(&maf, 1.0, 25.0, &mut weights);
                let weights: Array1<f64> = weights.into();
                // MU
                let mu = run_logistic_regression_b(&x, &y);
                // PI
                let pi: Array1<f64> = mu.map(|mu_i| mu_i * (1.0 - mu_i)).into();
                // X1
                let x1 = make_mat_x1(&x);
                let rho_vec = vec![0.1, 0.3, 0.5, 0.7, 0.9];

                println!("\trun Skat-O test");
                let mut calculator = SkatOCalculator::new(x1, z, weights, y, mu, pi, rho_vec);
                calculator.multiply_weight_to_z();
                calculator.multiply_p_to_z();
                // do this step before z multiply l (which is called in get_pmin_and_quantiles)
                let t_params = calculator.get_t_params();
                let (_t, quantiles) = calculator.get_pmin_and_quantiles();
                let t_pval = calculator.get_t_p_value(&quantiles, &t_params);

                println!("t_pval = {}", t_pval);
            }
        }
    }
}
