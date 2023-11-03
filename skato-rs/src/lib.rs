pub mod cffi;
pub mod preprocess;
mod stats;

use ndarray_linalg::{Cholesky, Eigh};
pub use preprocess::*;

pub use ndarray::{aview1, Array1, Array2, ArrayView1, Axis};
use statrs::distribution::{ChiSquared, Continuous, ContinuousCDF};

pub fn make_mat_x1(mat_x: &Array2<f64>) -> Array2<f64> {
    let mut v = vec![];
    let ncol = mat_x.ncols();
    let nrows = mat_x.nrows();
    for row in mat_x.rows() {
        v.push(1.0);
        for e in row {
            v.push(*e);
        }
    }
    Array2::from_shape_vec([nrows, ncol + 1], v).expect("fail to convert vec to array2")
}
pub struct TParams {
    mu_q: f64,
    var_q: f64,
    _ker_q: f64,
    lambda_arr: Array1<f64>,
    var_remain: f64,
    df: f64,
    tau: Vec<f64>,
}

pub struct SkatOCalculator {
    // covariable with additional columns for intercepts
    pub x1: Array2<f64>,
    // genotype data
    pub z: Array2<f64>,
    // weights on SNPs
    pub w: Array1<f64>,
    // response with 0.0 and 1.1
    pub y: Array1<f64>,
    // predicted probability
    pub mu: Array1<f64>,
    // inverse of the variance of y_i
    pub pi: Array1<f64>,
    // q_burden
    q_burden: f64,
    // q_stat
    q_skat: f64,
    rho_vec: Vec<f64>,
    q_rho_vec: Vec<f64>,
}

impl SkatOCalculator {
    pub fn new(
        x1: Array2<f64>,
        z: Array2<f64>,
        w: Array1<f64>,
        y: Array1<f64>,
        mu: Array1<f64>,
        pi: Array1<f64>,
        rho_vec: Vec<f64>,
    ) -> Self {
        // instance with tempory q values
        let mut cal = Self {
            x1,
            z,
            w,
            y,
            mu,
            pi,
            q_burden: 0.0,
            q_skat: 0.0,
            rho_vec,
            q_rho_vec: vec![],
        };
        let q_burden = cal._get_q_burden();
        let q_skat = cal._get_q_skat();
        // update q_burden and q_stat
        {
            cal.q_burden = q_burden;
            cal.q_skat = q_skat;
        }
        // update q_rho_vec
        {
            let mut v = vec![];
            for rho in cal.rho_vec.iter() {
                let rho = *rho;
                let mut q_rho = rho * q_burden + (1.0 - rho) * q_skat;
                q_rho /= 2.0; // not sure why but see SKAT_Optimal_Get_Q
                v.push(q_rho);
            }
            cal.q_rho_vec = v;
        }
        // return
        cal
    }
    /// see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415556/ (equation 4)
    fn _get_q_burden(&self) -> f64 {
        let base = self
            .y
            .iter()
            .zip(self.mu.iter())
            .zip(self.z.axis_iter(Axis(0)))
            .map(|((yi, mui), zi)| {
                let part1 = *yi - *mui;
                let part2 = self
                    .w
                    .iter()
                    .zip(zi.iter())
                    .map(|(wj, zij)| *wj * *zij)
                    .sum::<f64>();
                part1 * part2
            })
            .sum::<f64>();
        base * base
    }

    /// see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3415556/ (equation 6)
    fn _get_q_skat(&self) -> f64 {
        self.w
            .iter()
            .zip(self.z.axis_iter(Axis(1)))
            .map(|(wj, zj)| {
                let part1 = wj * wj;
                let mut part2 = zj
                    .iter()
                    .zip(self.y.iter())
                    .zip(self.mu.iter())
                    .map(|((zij, yi), mui)| *zij * (*yi - *mui))
                    .sum::<f64>();
                part2 = part2 * part2;
                part1 * part2
            })
            .sum::<f64>()
    }

    pub fn multiply_weight_to_z(&mut self) {
        self.z
            .columns_mut()
            .into_iter()
            .zip(self.w.iter())
            .for_each(|(mut zj, wj)| zj *= *wj)
    }

    /// NOTE: to improve computation efficiency, SKAT author did some
    /// similarity transformation:
    /// Given P1/2 * Kp * P1/2  is the taget matrix,
    /// let Q = P1/2 and Z0 = ZW (weighted Z)
    /// then Kp = Z * W * Rp * W * t(Z) = Z0 * Rp * t(Z0)
    /// Via cholesky decomposition Rp = L * L'
    /// the target matrix is now
    ///    Q * Z0 * Rp * t(Z0) * Q = Q * Z0 * L * L' * Z0' * Q
    /// which is a similar matrix to
    ///    Q * (Q * Z0 * L * L' * Z0' * Q ) * Q-1
    ///      = (Q * Q) * (Z0 * L * L' * Z0') * Q * Q-1
    ///      = P * (Z0 * L * L' * Z0')
    ///
    /// let M = diag (sqrt(pi)
    ///  NOTE: it can be prove that
    ///    P * M * M * P = P
    /// thus, P * (Z0 * L * L' * Z0')
    ///     = P * M * M * P * (Z0 * L * L' * Z0')
    ///     = M * P * Z0 * L * L' * Z0' * P * M
    ///
    /// This is more efficient as
    /// 1. We don't have to calculate P1/2 , and
    /// 2. We save the result the intermediate result t(Z0)*P*Z0 for different
    /// Rp

    /// Relevant R code:
    ///
    /// Z <- t(t(Z) * (weights))
    /// Z1 <- (Z * sqrt(pi_1)) - (X1 * sqrt(pi_1)) %*% solve(t(X1) %*% (X1 * pi_1)) %*%
    /// (t(X1) %*% (Z * pi_1))

    /// for each rho
    ///   r.corr <- r.all[i]
    ///   R.M <- diag(rep(1 - r.corr, p.m)) + matrix(rep(r.corr, p.m * p.m), ncol = p.m)
    ///   L <- chol(R.M, pivot = TRUE)
    ///   Z2 <- Z1 %*% t(L)
    ///   K1 <- t(Z2) %*% Z2
    ///   # eigenvalues
    ///   lambda.all[[i]] <- Get_Lambda(K1)

    pub fn multiply_p_to_z(&mut self) {
        use ndarray_linalg::solve::Inverse;
        // column vector
        let pi = self.pi.view();
        let z = self.z.view();
        let x1 = self.x1.view();

        // Z1 <- (Z * sqrt(pi_1)) -
        //        a
        // (X1 * sqrt(pi_1)) %*%
        //   b
        // solve(t(X1) %*% (X1 * pi_1)) %*%
        //  c
        // (t(X1) %*% (Z * pi_1))
        // d

        let mut a = z.to_owned();
        a.axis_iter_mut(Axis(0))
            .zip(pi.iter())
            .for_each(|(mut z_i, pi_i)| z_i *= pi_i.sqrt());
        let mut b = x1.to_owned();
        b.axis_iter_mut(Axis(0))
            .zip(pi.iter())
            .for_each(|(mut x_i, pi_i)| x_i *= pi_i.sqrt());
        let mut c = x1.to_owned();
        c.axis_iter_mut(Axis(0))
            .zip(pi.iter())
            .for_each(|(mut x_i, pi_i)| x_i *= *pi_i);
        let c = x1.t().dot(&c).inv().expect("fail to invert");
        let mut d = z.to_owned();
        d.axis_iter_mut(Axis(0))
            .zip(pi.iter())
            .for_each(|(mut z_i, pi_i)| z_i *= *pi_i);
        let d = x1.t().dot(&d);
        self.z = a - b.dot(&c).dot(&d);
    }

    // calculate z1 without changing self.z. i.e. z1 is calculated on the fly
    // for each rho
    pub fn multiply_l_to_z(&self, rho: f64) -> Array2<f64> {
        // prepare r_rho matrix
        let nsnp = self.w.len();
        let mut r_rho = Array2::from_elem([nsnp, nsnp], rho);
        r_rho.diag_mut().fill(1.0);
        // skat-o R code L is upper triangular matrix and then use tranpose
        // here we direct get the lower triangular matrix
        let l = r_rho.cholesky(ndarray_linalg::UPLO::Lower).unwrap();
        self.z.dot(&l).mapv_into(|x| x / 2.0f64.sqrt())
    }

    pub fn _get_sorted_nonneg_lambdas(mat: &Array2<f64>) -> Array1<f64> {
        //  out.s <- eigen(K, symmetric = TRUE, only.values = TRUE)
        //  lambda1 <- out.s$values
        //  IDX1 <- which(lambda1 >= 0)
        //  # eigenvalue bigger than sum(eigenvalues)/1000
        //  IDX2 <- which(lambda1 > mean(lambda1[IDX1]) / 1e+05)
        //  lambda <- lambda1[IDX2]
        let (eigvals, _eigvec) = mat.eigh(ndarray_linalg::UPLO::Lower).unwrap();
        let mut eigvals: Vec<_> = eigvals.into_iter().filter(|x| *x >= 0.0).collect();
        let thres = eigvals.iter().sum::<f64>() / eigvals.len() as f64 / 1e+5;
        eigvals.retain(|x| *x > thres);
        eigvals.sort_by(|a, b| a.partial_cmp(b).unwrap().reverse());
        eigvals.into()
    }

    pub fn _get_q_rho_params(&self, rho: f64) -> (f64, f64, f64) {
        // calculate z1 and lambdas for each q_skat
        let z1 = self.multiply_l_to_z(rho);

        let lambda = Self::_get_sorted_nonneg_lambdas(&z1.t().dot(&z1));
        let c1: f64 = lambda.sum();
        let c2: f64 = lambda.map(|x| *x * *x).sum();
        let c3: f64 = lambda.map(|x| *x * *x * *x).sum();
        let c4: f64 = lambda.map(|x| *x * *x * *x * *x).sum();

        // get mu/sigma/df
        // see Get_Liu_Params_Mod
        let mu_q = c1;
        let sigma_q = (2.0 * c2).sqrt();
        let df = {
            let sym = c3 / c2.powf(3.0 / 2.0);
            let kert = c4 / (c2 * c2);
            match sym * sym > kert {
                true => {
                    let a = 1.0 / (sym - (sym * sym - kert).sqrt());
                    let d = sym * a * a * a - a * a;
                    a * a - 2.0 * d
                }
                false => 1.0 / kert,
            }
        };
        (mu_q, sigma_q, df)
    }

    pub fn get_pmin_and_quantiles(&self) -> (f64, Vec<f64>) {
        let mut params = Vec::<(f64, f64, f64)>::with_capacity(self.q_rho_vec.len());
        let mut quantiles = Vec::<f64>::with_capacity(self.q_rho_vec.len());
        let mut pval_min = 1.0;
        //assert!(false);
        // calculate params and pval for each q_skat and get the min
        for (q_rho, rho) in self.q_rho_vec.iter().zip(self.rho_vec.iter()) {
            let (mu_q, sigma_q, df) = self._get_q_rho_params(*rho);
            params.push((mu_q, sigma_q, df));
            // Lee2012AJHG, equation 11
            let q_norm = (q_rho - mu_q) * (2.0 * df).sqrt() / sigma_q + df;
            let pval = 1.0 - ChiSquared::new(df).unwrap().cdf(q_norm);
            if pval < pval_min {
                pval_min = pval;
            }
        }
        // find quantile of 1-T in the chisq-mixture distribution for each q_stat
        let t = pval_min;
        for (mu_q, sigma_q, df) in params {
            let quantile_orig = ChiSquared::new(df).unwrap().inverse_cdf(1.0 - t);
            let quantile_q_stat = (quantile_orig - df) / (2.0 * df).sqrt() * sigma_q + mu_q;
            quantiles.push(quantile_q_stat);
        }
        (t, quantiles)
    }

    // get param for distribution for t
    // - mu_q
    // - var_q
    // - ker_q
    // - lambda
    // - var_remain
    // - df
    // - tau: dim=(5,)
    pub fn get_t_params(&self) -> TParams {
        let z1 = self.z.mapv(|zij| zij / 2.0f64.sqrt());
        use ndarray::prelude::*;
        // row mean
        let z_mean = z1.mean_axis(Axis(1)).unwrap();
        // cof1
        let zmean_sq_sum: f64 = z_mean.iter().map(|e| *e * *e).sum();
        let cof1: Array1<f64> = z1
            .columns()
            .into_iter()
            .map(|col| col.dot(&z_mean) / zmean_sq_sum)
            .collect();
        // item1

        let mut item1 = Array2::zeros([z1.nrows(), z1.ncols()]);
        item1
            .rows_mut()
            .into_iter()
            .zip(z_mean.iter())
            .for_each(|(mut item1_i, z_mean_i)| {
                item1_i.assign(&cof1);
                item1_i *= *z_mean_i;
            });

        let item2 = &z1 - &item1;

        let lambda_arr = Self::_get_sorted_nonneg_lambdas(&item2.t().dot(&item2));
        let var_remain = (item1.t().dot(&item1) * (&item2.t().dot(&item2))).sum() * 4.0;

        let mu_q: f64 = lambda_arr.iter().sum();
        let c2 = lambda_arr.iter().map(|e| *e * *e).sum::<f64>();
        let c4 = lambda_arr.iter().map(|e| e.powi(4)).sum::<f64>();
        let var_q = c2 * 2.0 + var_remain;
        let ker_q = c4 / c2 / c2 * 12.0;
        let df = 12.0 / ker_q;

        let nsnp = z1.ncols() as f64;
        let tau: Vec<f64> = self
            .rho_vec
            .iter()
            .map(|rho| {
                let term1 =
                    nsnp * nsnp * *rho + cof1.iter().map(|x| *x * *x).sum::<f64>() * (1.0 - *rho);
                term1 * zmean_sq_sum
            })
            .collect();

        TParams {
            mu_q,
            var_q,
            _ker_q: ker_q,
            lambda_arr,
            var_remain,
            df,
            tau,
        }
    }

    /// see SKAT R package function `SKAT_Optimal_Integrate_Func_Davies`
    pub fn _t_stat_pdf_davies(
        x: f64,
        quantiles: ArrayView1<f64>,
        rho_arr: ArrayView1<f64>,
        t_params: &TParams,
    ) -> f64 {
        let temp1 = &aview1(t_params.tau.as_slice()) * x;
        // TODO: need to improve
        let temp = (&quantiles - temp1) / (1.0 - &rho_arr);
        let mut temp_min = f64::MAX;
        temp.iter().for_each(|x| {
            if *x < temp_min {
                temp_min = *x;
            }
        });

        let temp = if temp_min > t_params.lambda_arr.sum() * 1e4 {
            0.0
        } else {
            let min1_temp = temp_min - t_params.mu_q;
            let sd1 = (t_params.var_q - t_params.var_remain).sqrt() / t_params.var_q.sqrt();
            let min1_st = min1_temp * sd1 + t_params.mu_q;

            // see SKAT_davies function
            // asserts
            let nlambda = t_params.lambda_arr.len();
            let mut trace_arr = [0.0f64; 7];
            let mut idx_fault = 0i32;
            let mut res = 0.0f64;
            cffi::qfc(
                t_params.lambda_arr.as_slice().unwrap(), // lb1
                &vec![0.0; nlambda],                     // nc1
                &vec![1i32; nlambda],                    // n1
                &(nlambda as i32),                       // r1
                &0.0f64,                                 // sigma,
                &min1_st,                                // c1
                &10000i32,                               // lim1
                &(1e-4),                                 // acc
                &mut trace_arr,                          // trace_arr
                &mut idx_fault,                          // idx_fault
                &mut res,                                // result
            )
            .unwrap();
            assert_eq!(idx_fault, 0);
            res = 1.0 - res;
            if res > 1.0 {
                res = 1.0;
            }
            res
        };

        let p = (1.0 - temp) * ChiSquared::new(1.0).unwrap().pdf(x);
        p
    }

    pub fn _t_stat_pdf_liu(
        x: f64,
        quantiles: ArrayView1<f64>,
        rho_arr: ArrayView1<f64>,
        t_params: &TParams,
    ) -> f64 {
        let temp1 = &aview1(t_params.tau.as_slice()) * x;
        // TODO: need to improve
        let temp = (&quantiles - temp1) / (1.0 - &rho_arr);
        let mut temp_min = f64::MAX;
        temp.iter().for_each(|x| {
            if *x < temp_min {
                temp_min = *x;
            }
        });
        let temp_q = (temp_min - t_params.mu_q) / t_params.var_q.sqrt()
            * (2.0 * t_params.df).sqrt()
            + t_params.df;

        let p = ChiSquared::new(t_params.df).unwrap().cdf(temp_q)
            * ChiSquared::new(1.0).unwrap().pdf(x);
        p
    }

    pub fn get_t_p_value(&self, quantiles: &Vec<f64>, t_params: &TParams) -> f64 {
        use peroxide::numerical::integral::{integrate, Integral};
        let rho_arr = aview1(self.rho_vec.as_slice());
        let quantiles = aview1(quantiles.as_slice());
        let re = integrate(
            // |x| Self::_t_stat_pdf_liu(x, quantiles, rho_arr, t_params),
            |x| Self::_t_stat_pdf_davies(x, quantiles, rho_arr, t_params),
            (0.0, 40.0),
            Integral::G7K15(1e-8, 30),
        );
        1.0 - re
    }
}

#[test]
/// Result should resemble results from the following R code
///     library(SKAT)
///     data("SKAT.example")
///     attach(SKAT.example)
///     obj <- SKAT_Null_Model(y.b ~ X, out_type = "D")
///     out.b <- SKAT(Z, obj, r.corr = c(0.1, 0.3, 0.5, 0.7, 0.9))
///     print(out.b$p.value)
fn test_skato_example() {
    use ndarray::arr2;
    // load dataset
    let (mat_x, mat_z, vec_yb) = load_test_datasets();
    // get MAF
    let maf = get_maf(&mat_z);
    // get weights
    let mut weights = Vec::new();
    get_beta_weights(&maf, 1.0, 25.0, &mut weights);
    let mu = run_logistic_regression_b(&mat_x, &vec_yb);
    // save_matrix(&mu.clone().into_shape([mu.len(), 1]).unwrap(), "out.mu.txt");
    let y: Array1<_> = vec_yb.into_raw_vec().into();

    // build calculator instance
    let x1 = make_mat_x1(&mat_x).into();
    let z = mat_z;
    let w = weights.into();
    let pi: Array1<f64> = mu.map(|mu_i| mu_i * (1.0 - mu_i)).into();
    let mu = mu.into();
    let rho_vec = vec![0.1, 0.3, 0.5, 0.7, 0.9];
    let mut calculator = SkatOCalculator::new(x1, z, w, y, mu, pi, rho_vec);
    calculator.multiply_weight_to_z();
    calculator.multiply_p_to_z();
    // do this step before z multiply l (which is called in get_pmin_and_quantiles)
    let t_params = calculator.get_t_params();
    let (_t, quantiles) = calculator.get_pmin_and_quantiles();
    let t_pval = calculator.get_t_p_value(&quantiles, &t_params);

    assert!((t_pval - 0.0843).abs() < 0.001);
    println!("t_pval = {}", t_pval);
}
