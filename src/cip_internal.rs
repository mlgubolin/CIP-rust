#[derive(Clone)]
pub struct GridParameters {
    f: Vec<f64>,
    u: Vec<f64>,
    g: Vec<f64>,
    df_star: Vec<f64>,
    dx: f64,
    dt: f64,
}

impl GridParameters {
    pub fn new(f: Vec<f64>, u: Vec<f64>, g: Vec<f64>, dx: f64, dt: f64) -> GridParameters {
        return GridParameters {
            f: f,
            u: u,
            g: g,
            dx: dx,
            dt: dt,
            df_star: vec![],
        };
    }
}

#[derive(Clone)]
pub struct Grids {
    pub previous_previous: GridParameters,
    pub previous: GridParameters,
    pub current: GridParameters,
    pub next: GridParameters,
    pub next_next: GridParameters,
}

impl Grids {
    pub fn set_grids(
        previous_previous: GridParameters,
        previous: GridParameters,
        current: GridParameters,
        next: GridParameters,
        next_next: GridParameters,
    ) -> Grids {
        return Grids {
            previous_previous: previous_previous,
            previous: previous,
            current: current,
            next: next,
            next_next: next_next,
        };
    }
}

struct FStar {
    previous_previous: Vec<f64>,
    previous: Vec<f64>,
    current: Vec<f64>,
    next: Vec<f64>,
}

impl FStar {
    fn new(
        previous_previous: Vec<f64>,
        previous: Vec<f64>,
        current: Vec<f64>,
        next: Vec<f64>,
    ) -> FStar {
        return FStar {
            previous_previous: previous_previous,
            previous: previous,
            current: current,
            next: next,
        };
    }
}

pub fn cip_calculation(grids: Grids) -> Vec<f64> {
    let f_star = calculate_f_star(&grids);
    let df_star = calculate_df_star(&grids, &f_star);

    let h_vector = h(&f_star);
    let g_vector = g(&df_star.next, &df_star.previous, &grids.current.dx);

    let df_future_next = estimate_df_plus(
        &grids.next,
        &f_star.current,
        &f_star.previous,
        &df_star.current,
        &df_star.previous,
    );

    let df_future_previous = estimate_df_plus(
        &grids.previous,
        &f_star.previous_previous,
        &f_star.previous,
        &grids.previous_previous.df_star,
        &df_star.previous,
    );

    let g_plus_vector = g(&df_future_next, &df_future_previous, &grids.current.dx);

    let delta_f_plus = delta_f(
        &grids,
        &df_star.next,
        &df_star.current,
        &f_star.next,
        &f_star.current,
    );

    let delta_f_minus = delta_f(
        &grids,
        &df_star.current,
        &df_star.previous,
        &f_star.current,
        &f_star.previous,
    );

    return calculate_final_solution(
        h_vector,
        g_vector,
        g_plus_vector,
        delta_f_plus,
        delta_f_minus,
        df_future_next,
        df_future_previous,
        grids.current.dx,
    );
}

pub fn cip_initialization(mut first_grid_set: Grids) {
    let f_star = calculate_f_star(&first_grid_set);
    let df_star = calculate_df_star(&first_grid_set, &f_star);
    let df_previous_previous = vec![0.0; first_grid_set.current.f.len()];

    first_grid_set.previous.df_star = estimate_df_plus(
        &first_grid_set.previous,
        &f_star.previous_previous,
        &f_star.previous,
        &df_previous_previous,
        &df_star.previous,
    );
}

fn calculate_f_star(grids: &Grids) -> FStar {
    let vector_lenght = grids.current.f.len();
    let mut f_star_previous_previous = Vec::new();
    let mut f_star_previous = Vec::new();
    let mut f_star_current = Vec::new();
    let mut f_star_next = Vec::new();

    for i in 0..vector_lenght {
        f_star_previous_previous.push(
            grids.previous_previous.f[i]
                + grids.previous_previous.g[i] * grids.previous_previous.dt,
        );
        f_star_previous.push(grids.previous.f[i] + grids.previous.g[i] * grids.previous.dt);
        f_star_current.push(grids.current.f[i] + grids.current.g[i] * grids.current.dt);
        f_star_next.push(grids.next.f[i] + grids.next.g[i] * grids.next.dt);
    }

    return FStar::new(
        f_star_previous_previous,
        f_star_previous,
        f_star_current,
        f_star_next,
    );
}

fn calculate_df_star(grids: &Grids, f_star: &FStar) -> FStar {
    let df_star_previous_previous = &grids.previous_previous.df_star;

    let df_star_previous = calculate_df_star_individual_vector(
        &grids.previous_previous,
        &grids.previous,
        &grids.current,
        &f_star,
    );

    let df_star_current =
        calculate_df_star_individual_vector(&grids.previous, &grids.current, &grids.next, &f_star);

    let df_star_next =
        calculate_df_star_individual_vector(&grids.current, &grids.next, &grids.next_next, &f_star);

    return FStar::new(
        df_star_previous_previous.to_vec(),
        df_star_previous,
        df_star_current,
        df_star_next,
    );
}
fn calculate_df_star_individual_vector(
    previous: &GridParameters,
    current: &GridParameters,
    next: &GridParameters,
    f_star: &FStar,
) -> Vec<f64> {
    let vector_lenght = current.f.len();
    let f_derivative_estimate = estimate_derivative(&previous, &current, &next);

    let mut df_star = Vec::new();

    for i in 0..vector_lenght {
        df_star.push(
            f_derivative_estimate[i]
                + (f_star.next[i] - f_star.previous[i] - next.f[i] + previous.f[i])
                    / (2.0 * current.dx),
        );
    }

    return df_star;
}

fn estimate_derivative(
    previous: &GridParameters,
    current: &GridParameters,
    next: &GridParameters,
) -> Vec<f64> {
    let vector_lenght = current.f.len();
    let mut f_derivative_estimate = Vec::new();

    for i in 0..vector_lenght {
        let next_deriv = (next.f[i] - current.f[i]) / current.dx;
        let prev_deriv = (current.f[i] - previous.f[i]) / current.dx;

        f_derivative_estimate.push((next_deriv + prev_deriv) / 2.0);
    }

    return f_derivative_estimate;
}

fn h(f_star: &FStar) -> Vec<f64> {
    let vector_lenght = f_star.current.len();
    let mut h_vector = Vec::new();

    for i in 0..vector_lenght {
        h_vector.push(
            (18.0 * f_star.next[i] + 156.0 * f_star.current[i] + 18.0 * f_star.previous[i]) / 192.0,
        );
    }
    return h_vector;
}

// fn H_mean(
//     current: &GridParameters,
//     f_star: &FStar,
//     // }
//     df_star: &FStar,
//     dx: &f64
// )-> Vec<f64>{
//     let vector_lenght = f_star.current.len();
//     let mut h_mean_vector = Vec::new();

//     let (a_vec, b_vec, x_vec) = calculate_coefficients(current, f_star, df_star);

//     for i in 0..vector_lenght{
//         h_mean_vector.push(
//             3.0 * a_vec[i] * x_vec[i].powi(2)
//           + 2.0 * b_vec[i] * x_vec[i]
//           + df_star.previous[i]
//         );
//     }
//    return h_mean_vector;
// }

fn estimate_df_plus(
    current: &GridParameters,
    f_star_previous: &Vec<f64>,
    f_star_current: &Vec<f64>,
    df_star_previous: &Vec<f64>,
    df_star_current: &Vec<f64>,
) -> Vec<f64> {
    let vector_lenght = current.f.len();
    let mut df_estimate = Vec::new();

    let (x_vec, a_vec, b_vec) = calculate_coefficients(
        &f_star_previous,
        &f_star_current,
        &df_star_previous,
        &df_star_current,
        current,
    );

    for i in 0..vector_lenght {
        df_estimate.push(
            3.0 * a_vec[i] * x_vec[i].powi(2) + 2.0 * b_vec[i] * x_vec[i] + df_star_previous[i],
        );
    }
    return df_estimate;
}

fn calculate_coefficients(
    f_star_previous: &Vec<f64>,
    f_star_current: &Vec<f64>,
    df_star_previous: &Vec<f64>,
    df_star_current: &Vec<f64>,
    current: &GridParameters,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let vector_lenght = current.f.len();
    let dx = current.dx;
    let dt = current.dt;
    let mut x_vec = Vec::new();
    let mut a_vec = Vec::new();
    let mut b_vec = Vec::new();

    for i in 0..vector_lenght {
        x_vec.push(dx - current.u[i] * dt);

        a_vec.push(
            (df_star_current[i] * df_star_previous[i]) / dx.powi(2)
                - 2.0 * (f_star_current[i] - f_star_previous[i]) / dx.powi(3),
        );

        b_vec.push(
            3.0 * (f_star_current[i] - f_star_previous[i]) / dx.powi(2)
                - (df_star_current[i] + 2.0 * df_star_previous[i]) / dx,
        );
    }

    return (x_vec, a_vec, b_vec);
}

fn g(next: &Vec<f64>, previous: &Vec<f64>, dx: &f64) -> Vec<f64> {
    let vector_lenght = next.len();
    let mut g_vector = Vec::new();

    for i in 0..vector_lenght {
        g_vector.push(5.0 * (previous[i] - next[i]) / (192.0 * dx));
    }
    return g_vector;
}

fn delta_f(
    grids: &Grids,
    df_next: &Vec<f64>,
    df_current: &Vec<f64>,
    f_next: &Vec<f64>,
    f_current: &Vec<f64>,
) -> Vec<f64> {
    let vector_lenght = f_current.len();
    let mut delta_f = Vec::new();
    let dx = grids.current.dx;

    for i in 0..vector_lenght {
        let kappa = grids.current.dt * (grids.next.u[i] - grids.current.u[i]) / (2.0 * dx.powi(2)); //Need to check the validity of this riemman problem
        delta_f.push(
            (-kappa / 8.0 + kappa.powi(2) / 8.0 + kappa.powi(3) / 6.0 - kappa.powi(4) / 4.0)
                * df_next[i]
                * dx.powi(2)
                + (kappa / 8.0 + kappa.powi(2) / 8.0 - kappa.powi(3) / 6.0 - kappa.powi(4) / 4.0)
                    * df_current[i]
                    * dx.powi(2)
                + (kappa / 2.0 - 3.0 * kappa.powi(2) / 4.0 - kappa.powi(4) / 4.0) * f_next[i] * dx
                + (kappa / 2.0 - 3.0 * kappa.powi(2) / 4.0 - kappa.powi(4) / 4.0)
                    * f_current[i]
                    * dx,
        );
    }

    return delta_f;
}

fn calculate_final_solution(
    h_vector: Vec<f64>,
    g_vector: Vec<f64>,
    g_plus_vector: Vec<f64>,
    delta_f_plus: Vec<f64>,
    delta_f_minus: Vec<f64>,
    df_future_next: Vec<f64>,
    df_future_previous: Vec<f64>,
    dx: f64,
) -> Vec<f64> {
    let vector_lenght = h_vector.len();
    let mut final_solution = Vec::new();

    for i in 0..vector_lenght {
        let right_sum = h_vector[i] + g_vector[i] - (delta_f_plus[i] - delta_f_minus[i]) / dx;
        final_solution.push(
            (192.0 * (right_sum - g_plus_vector[i])
                - 18.0 * df_future_next[i]
                - 18.0 * df_future_previous[i])
                / 156.0,
        )
    }

    return final_solution;
}
