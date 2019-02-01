// This crate is a library
#![crate_type = "lib"]
// The library is named "rary"
#![crate_name = "cip_rust"]
mod cip_internal;


// Returns f^{n+1}_i
pub fn CIP_1D(grids: cip_internal::Grids) -> Vec<f64>{
    return cip_internal::calculation(grids);
}