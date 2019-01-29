mod cip_internal;
use std::vec::Vec;

pub struct grid_setup {
    current_grid: grid_variables,
    previous_grid: grid_variables,
    next_grid: grid_variables,
}

struct grid_variables {
    f: Vec<f64>,
    u: Vec<f64>,
    g: Vec<f64>,
    dx: f64,
    dt: f64,
}

// Returns f^{n+1}_i
pub fn cip_1d(grid_setup: grid_setup) -> Vec<f64> {
    return vec![0.1, 1.1];
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
