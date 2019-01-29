mod cip_internal;
use std::vec::Vec;

pub struct Grids {
    current_grid: GridParameters,
    previous_grid: GridParameters,
    next_grid: GridParameters,
}

impl Grids {
    fn set_current_grid(f: Vec<f64>, u: Vec<f64>, g: Vec<f64>, dx: f64, dt: f64) -> GridParameters {
        return GridParameters {
            f: f,
            u: u,
            g: g,
            dx: dx,
            dt: dt,
        };
    }
}

struct GridParameters {
    f: Vec<f64>,
    u: Vec<f64>,
    g: Vec<f64>,
    dx: f64,
    dt: f64,
}

// Returns f^{n+1}_i
pub fn cip_1d(Grids: Grids) -> Vec<f64> {
    return vec![0.1, 1.1];
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
