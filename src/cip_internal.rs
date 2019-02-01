pub struct GridParameters {
    f: Vec<f64>,
    u: Vec<f64>,
    g: Vec<f64>,
    df: Vec<f64>,
    dx: f64,
    dt: f64,
}

impl GridParameters{
    pub fn new(
        f: Vec<f64>,
        u: Vec<f64>,
        g: Vec<f64>,
        dx: f64,
        dt: f64)
         -> GridParameters{
        
        return GridParameters{
            f: f,
            u: u,
            g: g,
            dx: dx,
            dt: dt,
            df: vec![]
        }
    }
}


pub struct Grids {
    previous_previous: GridParameters,
    previous: GridParameters,
    current: GridParameters,
    next: GridParameters,
    next_next: GridParameters,
}

impl Grids {
    pub fn set_grid_parameters(f: Vec<f64>, u: Vec<f64>, g: Vec<f64>, dx: f64, dt: f64) -> GridParameters {
        return GridParameters::new(f, u, g, dx, dt);
    }
}


struct FStar{
    previous: Vec<f64>,
    current: Vec<f64>,
    next: Vec<f64>,
}

impl FStar{
    fn new(previous: Vec<f64>,
    current: Vec<f64>,
    next: Vec<f64>) -> FStar{
        return FStar{
            previous: previous,
            current: current,
            next: next,
        }
    }

}

pub fn calculation(grids: Grids) -> Vec<f64>{

    let f_star = calculate_f_star(&grids);
    let df_star = calculate_df_star(&grids, &f_star);

    let h_vector = H(&f_star);
    let h_mean_vector = H_mean(&f_star);
    let g_vector = G(
        &f_star.next,
        &f_star.previous,
        &grids.current.dx
    );

    

    let delta_F_plus = deltaF(
        &grids,
        &df_star.next, 
        &df_star.current,
        &f_star.next,
        &f_star.current
    );

    let delta_F_minus = deltaF(
        &grids,
        &df_star.current, 
        &df_star.previous,
        &f_star.current,
        &f_star.previous
    );

    return calculate_next_f(
        h_vector,
        h_mean_vector,
        g_vector,
        g_plus_vector,
        delta_F_plus,
        delta_F_minus,
        grids.current.dx
    )
}

fn calculate_f_star(grids: &Grids) -> FStar{
    let vector_lenght = grids.current.f.len();
    let mut f_star_previous = Vec::new();
    let mut f_star_current = Vec::new();
    let mut f_star_next = Vec::new();

    for i in 0..vector_lenght{
       f_star_previous.push(grids.previous.f[i] + grids.previous.g[i] * grids.previous.dt);
       f_star_current.push(grids.current.f[i] + grids.current.g[i] * grids.current.dt);
       f_star_next.push(grids.next.f[i] + grids.next.g[i] * grids.next.dt);
    }

    return FStar::new(f_star_previous, f_star_current, f_star_next);
}

fn calculate_df_star(grids: &Grids,f_star: &FStar) -> FStar{
    let df_star_previous = calculate_df_star_individual_vector(
           &grids.previous_previous, 
           &grids.previous, 
           &grids.current, 
           &f_star);
 
    let df_star_current = calculate_df_star_individual_vector(
           &grids.previous, 
           &grids.current, 
           &grids.next, 
           &f_star);

    let df_star_next = calculate_df_star_individual_vector(
           &grids.current, 
           &grids.next, 
           &grids.next_next, 
           &f_star);

    return FStar::new(
        df_star_previous,
        df_star_current,
        df_star_next
    );
}

fn calculate_df_star_individual_vector(previous: &GridParameters, 
  current: &GridParameters,
  next: &GridParameters,
  f_star: &FStar)
   -> Vec<f64>{
    let vector_lenght = current.f.len();
    let f_derivative_estimate = estimate_derivative(&previous, &current, &next);
    let mut df_star = Vec::new();

    for i in 0..vector_lenght{
       df_star.push(f_derivative_estimate[i] 
        + (f_star.next[i] - f_star.previous[i] - next.f[i] + previous.f[i])
        / (2.0 * current.dx));
    }

    return df_star;
}

fn estimate_derivative(previous: &GridParameters, 
  current: &GridParameters,
  next: &GridParameters) -> Vec<f64>{
    let vector_lenght = current.f.len();
    let mut f_derivative_estimate = Vec::new();

    for i in 0..vector_lenght{
        let next_deriv = (next.f[i] - current.f[i])/current.dx;
        let prev_deriv = (current.f[i] - previous.f[i])/current.dx;

        f_derivative_estimate.push((next_deriv + prev_deriv) / 2.0);
    }

    return f_derivative_estimate;
}

    // let vector_lenght = grids.current.f.len()
    // let (x_vec, a_vec, b_vec) = calculate_coefficients(&grids.current, f_star);

    // for i in 0..vector_lenght{
    //    df_star.push(3.0 * a_vec[i] * x_vec[i].powi(2)
    //     + 2.0 * b_vec[i] * x_vec[i] 
    //     + previous.df[i]);
    // }

fn H(f_star: &FStar)-> Vec<f64>{
    let vector_lenght = f_star.current.len();
    let mut h_vector = Vec::new();

    for i in 0..vector_lenght{
        h_vector.push(
            (18.0 * f_star.next[i]
          + 156.0 * f_star.current[i]
          + 18.0 * f_star.previous[i])
          /192.0);
    }
   return h_vector;
}

fn H_mean(
    current: &GridParameters,
    f_star: &FStar,
    df_star: &FStar,
    dx: &f64
)-> Vec<f64>{
    let vector_lenght = f_star.current.len();
    let mut h_mean_vector = Vec::new();

    let (a_vec, b_vec, x_vec) = calculate_coefficients(current, f_star, df_star);

    for i in 0..vector_lenght{
        h_mean_vector.push(
            3.0 * a_vec[i] * x_vec[i].powi(2)
          + 2.0 * b_vec[i] * x_vec[i]
          + df_star.previous[i]
        );
    }
   return h_mean_vector;
}
fn estimate_f_plus(
    x_vec: &Vec<f64>,
    a_vec: &Vec<f64>,
    b_vec: &Vec<f64>,
    f_star: &FStar,
    dx: &f64
)-> Vec<f64>{
    let vector_lenght = f_star.current.len();
    let mut h_mean_vector = Vec::new();

    for i in 0..vector_lenght{
        h_mean_vector.push(
            3.0 * (f_star.next[i]
          + 2.0 * f_star.current[i]
          +       f_star.previous[i])
          /32.0);
    }
   return h_mean_vector;
}

fn calculate_coefficients(current: &GridParameters, f_star: &FStar, df_star: &FStar)-> (Vec<f64>, Vec<f64>, Vec<f64>){
    let vector_lenght = current.f.len();
    let mut x_vec = Vec::new();
    let mut a_vec = Vec::new();
    let mut b_vec = Vec::new();

    for i in 0..vector_lenght{
        x_vec.push(current.dx 
            - current.u[i] * current.dt);

        a_vec.push(
            (df_star.current[i] * df_star.previous[i])/current.dx.powi(2)
            - 2.0*(f_star.current[i] - f_star.previous[i])/current.dx.powi(3)
            );

        b_vec.push(3.0 * (f_star.current[i] - f_star.previous[i])/current.dx.powi(2)
            - (df_star.current[i] + 2.0 * df_star.previous[i])/current.dx
            );
    }

    return (x_vec, a_vec, b_vec);
}

fn G(next: &Vec<f64>, previous: &Vec<f64>, dx: &f64)-> Vec<f64>{
    let vector_lenght = next.len();
    let mut g_vector = Vec::new();

    for i in 0..vector_lenght{
        g_vector.push(
            5.0 * (previous[i]
            - next[i]) 
            / (192.0 * dx));
    }
   return g_vector;
}

fn deltaF(grids: &Grids, df_next: &Vec<f64>, df_current: &Vec<f64>,
  f_next: &Vec<f64>, f_current: &Vec<f64>)
  -> Vec<f64>{
    let vector_lenght = f_current.len();
    let mut delta_F = Vec::new();
    let dx = grids.current.dx;

    for i in 0..vector_lenght{
        let kappa = grids.current.dt
        * (grids.next.u[i] - grids.current.u[i])
        / (2.0 * dx.powi(2)); //Need to check the validity of this riemman problem
        delta_F.push(
            (-kappa / 8.0 + kappa.powi(2) / 8.0 + kappa.powi(3) / 6.0 - kappa.powi(4) / 4.0) * df_next[i] * dx.powi(2)
          + (kappa / 8.0 + kappa.powi(2) / 8.0 - kappa.powi(3) / 6.0 - kappa.powi(4) / 4.0) * df_current[i] * dx.powi(2)
          + (kappa/ 2.0 - 3.0 * kappa.powi(2) / 4.0 - kappa.powi(4) / 4.0) * f_next[i] * dx
          + (kappa/ 2.0 - 3.0 * kappa.powi(2) / 4.0 - kappa.powi(4) / 4.0) * f_current[i] * dx
        );
    }

    return delta_F;
}

fn calculate_next_f(
    h_vector: Vec<f64>,
    h_mean_vector: Vec<f64>,
    g_vector: Vec<f64>,
    g_plus_vector: Vec<f64>,
    delta_F_plus: Vec<f64>,
    delta_F_minus: Vec<f64>,
    dx: f64
) -> Vec<f64>{
    let vector_lenght = h_vector.len();
    let mut f_next = Vec::new();

    for i in 0..vector_lenght{
        f_next.push(
            h_vector[i]
          - h_mean_vector[i]
          + g_vector[i]
          - g_plus_vector[i]
          (delta_F_plus[i] - delta_F_minus[i])/dx
        )
    }
    return f_next;
}