use super::super::Commands;
use ishare::{genotype::rare::GenotypeRecords, indiv::Individuals, site::Sites};

pub fn main_skato(args: &Commands) {
    if let Commands::Skato { rec, genomes: _ } = args {
        let mut records = GenotypeRecords::from_parquet_file(rec);
        let ind_file = rec.with_extension("ind");
        let inds = Individuals::from_parquet_file(&ind_file);
        let sites_file = rec.with_extension("sit");
        let sites = Sites::from_parquet_file(&sites_file);
        println!("number of individuals: {}", inds.v().len());
        println!("number of sites      : {}", sites.len());
        // let freq_map = calc_allele_frequency(&mut records, inds.v().len() * 2, sites.len());

        // for (i, (k, v)) in freq_map.iter().enumerate() {
        //     if i < 5 {
        //         println!("k = {k}, v= {v}");
        //     }
        // }

        records.sort_by_position();
        let sitepos = sites.get_gw_pos_slice();
        println!(
            "site min_pos = {}, site max_pos = {}",
            sitepos[0],
            sitepos.last().unwrap()
        );
        let min_pos = sites.get_gw_pos_slice()[0];
        let max_pos = min_pos + 1000000;
        // first subset the records
        let all_records = records.records();
        let min_idx = all_records.partition_point(|x| x.get_position() < min_pos);
        let max_idx = all_records.partition_point(|x| x.get_position() < max_pos);
        let site_min_idx = sites.get_gw_pos_slice().partition_point(|x| *x < min_pos);
        let site_max_idx = sites.get_gw_pos_slice().partition_point(|x| *x < max_pos);
        let sel_records = &all_records[min_idx..max_idx];
        println!("number of selected records: {}", sel_records.len());
        println!(
            "number of selected sites  : {}",
            site_max_idx - site_min_idx
        );

        let nsites = site_max_idx - site_min_idx;
        let nsams = inds.v().len();

        {
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

            let mut calculator = SkatOCalculator::new(x1, z, weights, y, mu, pi, rho_vec);
            calculator.multiply_weight_to_z();
            calculator.multiply_p_to_z();
            // do this step before z multiply l (which is called in get_pmin_and_quantiles)
            let t_params = calculator.get_t_params();
            let (_t, quantiles) = calculator.get_pmin_and_quantiles();
            let t_pval = calculator.get_t_p_value(&quantiles, &t_params);

            println!("t_pval = {}", t_pval);
        }

        // run skato

        // after this we can focus on the full dataset
    }
}
