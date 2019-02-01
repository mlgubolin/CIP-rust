// This crate is a library
#![crate_type = "lib"]
// The library is named "rary"
#![crate_name = "cip_rust"]

use std::vec::Vec;
mod cip_internal;

// Returns f^{n+1}_i
pub fn cip_1d(grids: cip_internal::Grids) -> Vec<f64> {
    let FMinus = 1;
    let HPlus = cip_internal::calculate_GPlus(grids);
    return vec![0.1, 1.1];
}