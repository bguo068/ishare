extern crate libc;
use libc::{c_double, c_int};

extern "C" {
    fn qfc_1(
        lb1: *const c_double,
        nc1: *const c_double,
        n1: *const c_int,
        r1: *const c_int,
        sigma: *const c_double,
        c1: *const c_double,
        lim1: *const c_int,
        acc: *const c_double,
        trace_arr: *mut c_double,
        idx_fault: *mut c_int,
        res: *mut c_double,
    );
}

pub fn qfc(
    lb1: &[f64],
    nc1: &[f64],
    n1: &[i32],
    r1: &i32,
    sigma: &f64,
    c1: &f64,
    lim1: &i32,
    acc: &f64,
    trace_arr: &mut [f64],
    idx_fault: &mut i32,
    res: &mut f64,
) -> Result<(), &'static str> {
    // Safety checks can go here
    if lb1.len() == 0 || nc1.len() == 0 || trace_arr.len() == 0 {
        return Err("Input slices cannot be empty");
    }

    unsafe {
        qfc_1(
            lb1.as_ptr(),
            nc1.as_ptr(),
            n1.as_ptr(),
            r1,
            sigma,
            c1,
            lim1,
            acc,
            trace_arr.as_mut_ptr(),
            idx_fault,
            res,
        );
    }

    Ok(())
}

#[test]
fn test_qfc() {
    let lb1 = vec![0.0; 5];
    let nc1 = vec![0.0; 5];
    let n1 = vec![1i32, 5];
    let r1 = 5i32;
    let sigma = 0.0f64;
    let c1 = 0.0f64;
    let lim1 = 0i32;
    let acc = 0.0f64;
    let mut trace_arr = vec![0.0; 7];
    let mut idx_fault = 0i32;
    let mut res = 0.0f64;

    match qfc(
        &lb1,
        &nc1,
        &n1,
        &r1,
        &sigma,
        &c1,
        &lim1,
        &acc,
        &mut trace_arr,
        &mut idx_fault,
        &mut res,
    ) {
        Ok(()) => println!("Successfully called qfc_1"),
        Err(e) => println!("Error: {}", e),
    }
}
