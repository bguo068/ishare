use ndarray::{Array1, Array2, Axis};
use ndarray_glm::{Logistic, ModelBuilder};
use statrs::distribution::{Beta, Continuous};

use std::f64;
use std::fmt::{Debug, Display};
use std::str::FromStr;

pub fn load_matrix<'a, T: FromStr>(p: &str) -> Array2<T>
where
    <T as FromStr>::Err: Debug,
{
    let s = std::fs::read_to_string(p).expect("cannot read file");
    let mut ncol_line1 = 0;
    let mut nrow = 0;
    let mut a = vec![];

    for (i, line) in s.trim().split("\n").enumerate() {
        // number of columns for current line
        let mut ncol0 = 0;
        for col in line.split("\t") {
            ncol0 += 1;
            let num = col.parse::<T>().expect("can not parse fields");
            // add number to a vector
            a.push(num);
        }
        if i == 0 {
            ncol_line1 = ncol0;
        } else {
            // ensure number of columns is consistent between lines
            assert_eq!(ncol_line1, ncol0, "col number inconsistent between lines");
        }
        nrow += 1;
    }
    // non-copy conversion from Vec to Array2
    Array2::from_shape_vec([nrow, ncol_line1], a).expect("fail to convert Vec into Array2")
}

#[allow(dead_code)]
pub fn save_matrix<T: Display>(arr: &Array2<T>, p: &str) {
    use std::io::Write;
    let mut f = std::fs::File::create(p)
        .map(std::io::BufWriter::new)
        .unwrap();
    for row in arr.rows() {
        let mut count = 0;
        for i in row.iter() {
            if count == 0 {
                write!(f, "{}", *i).unwrap();
            } else {
                write!(f, "\t{}", *i).unwrap();
            }
            count += 1;
        }
        write!(f, "\n").unwrap();
    }
}

pub fn load_test_datasets() -> (Array2<f64>, Array2<f64>, Array1<f64>) {
    let mat_x = load_matrix::<f64>("./testdata/ex_data_x");
    let mat_z = load_matrix("./testdata/ex_data_z");
    // convert from Array2 to Array1
    let mat_yb = load_matrix("./testdata/ex_data_yb").into_raw_vec().into();
    (mat_x, mat_z, mat_yb)
}

pub fn get_maf(mat_z: &Array2<f64>) -> Vec<f64> {
    mat_z
        .axis_iter(Axis(1))
        .map(|col| {
            let mut ac = col.sum();
            let ta = (col.len() * 2) as f64;
            if ac * 2.0 > ta {
                ac = ta - ac;
            }
            ac / ta
        })
        .collect()
}

/// Get weights from a function of MAF. The function is the probability density
/// function of a beta distribution parametered by shape parameters a and b
pub fn get_beta_weights(maf: &[f64], a: f64, b: f64, weights: &mut Vec<f64>) {
    let d = Beta::new(a, b).expect("invalid beta parameters");
    weights.clear();
    let it = maf.iter().map(|maf| d.pdf(*maf));
    weights.extend(it);
}

/// run logistic regression using `ndarray-glm` crate.
/// Note: tried linfa-logistic v0.6.1 but it generated results very different
/// from python sklearn and R glm results
pub fn run_logistic_regression_b(mat_x: &Array2<f64>, vec_yb: &Array1<f64>) -> Array1<f64> {
    let y = &vec_yb;
    let model = ModelBuilder::<Logistic>::data(y, mat_x)
        .build()
        .expect("error in building model");
    let fit = model
        .fit_options()
        .tol(1e-6)
        .max_iter(100000)
        .fit()
        .expect("error when fitting model");
    let input = ndarray_glm::utility::one_pad(mat_x.view());
    let pred = fit.predict(&input, None);
    pred
}
