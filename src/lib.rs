// This crate is a library
#![crate_type = "lib"]
// The library is named "rary"
#![crate_name = "cip_rust"]
mod cip_internal;
pub type Grids = cip_internal::Grids;
pub type GridParameters = cip_internal::GridParameters;

// Returns f^{n+1}_i
pub fn cip_1d(grids: Grids) -> Vec<f64> {
    return cip_internal::cip_calculation(grids);
}

pub fn initialize(first_grid_set: Grids){
    return cip_internal::cip_initialization(first_grid_set);
}
