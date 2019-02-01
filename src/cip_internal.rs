#[derive(Clone)]
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

#[derive(Clone)]
pub struct Grids {
    current_grid: GridParameters,
    previous_grid: GridParameters,
    next_grid: GridParameters,
}

impl Grids {
    pub fn set_grid_parameters(f: Vec<f64>, u: Vec<f64>, g: Vec<f64>, dx: f64, dt: f64) -> GridParameters {
        return GridParameters::new(f, u, g, dx, dt);
    }
}

#[derive(Clone)]
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

pub fn calculation(grids: Grids){

    let f_star = calculate_fstar(&grids);
    let df_star = calculate_fplus(&grids, f_star.clone());

    let h_vector = H(&grids, df_star.clone());
    let g_vector = G(&grids, df_star.clone());


}

fn calculate_fstar(grids: &Grids) -> FStar{
    let vector_lenght = grids.current_grid.f.len();
    let mut f_star_previous = Vec::new();
    let mut f_star_current = Vec::new();
    let mut f_star_next = Vec::new();

    for i in 0..vector_lenght{
       f_star_previous.push(grids.previous_grid.f[i] + grids.previous_grid.g[i] * grids.previous_grid.dt as f64);
       f_star_current.push(grids.current_grid.f[i] + grids.current_grid.g[i] * grids.current_grid.dt  as f64);
       f_star_next.push(grids.next_grid.f[i] + grids.next_grid.g[i] * grids.next_grid.dt  as f64);
    }

    return FStar::new(f_star_previous, f_star_current, f_star_next);
}

fn calculate_fplus(grids: &Grids,f_star: FStar) -> Vec<f64>{
    let vector_lenght = grids.current_grid.f.len();
    let mut df_star = Vec::new();

    let (x_vec, a_vec, b_vec) = calculate_coefficients(grids.current_grid, f_star);

    for i in 0..vector_lenght{
       df_star.push(3.0 * a_vec[i] * x_vec[i].powi(2)
        + 2.0 * b_vec[i] * x_vec[i] 
        + grids.previous_grid.df[i]);
    }

    return df_star;

}

fn calculate_coefficients(current_grid: GridParameters, f_star: FStar)-> (Vec<f64>, Vec<f64>, Vec<f64>){
    let vector_lenght = current_grid.f.len();
    let mut x_vec = Vec::new();
    let mut a_vec = Vec::new();
    let mut b_vec = Vec::new();

    for i in 0..vector_lenght{
        x_vec.push(current_grid.dx 
            - current_grid.u[i] * current_grid.dt);

        a_vec.push(
            (f_star.current[i] * f_star.previous[i])/current_grid.dx.powi(2)
            - 2.0*(f_star.current[i] - f_star.previous[i])/current_grid.dx.powi(3)
            );

        b_vec.push(3.0 * (f_star.current[i] - f_star.previous[i])/current_grid.dx.powi(2)
            - (f_star.current[i] + 2.0 * f_star.previous[i])/current_grid.dx
            );
    }

    return (x_vec, a_vec, b_vec);
}

fn H(grids: &Grids, df_star: Vec<f64>)-> Vec<f64>{
    let vector_lenght = grids.current_grid.f.len();
    let mut h_vector = Vec::new();

    for i in 0..vector_lenght{
        h_vector.push((18.0*grids.next_grid.f[i] + 156.0*grids.current_grid.f[i] + 18.0*grids.previous_grid.f[i])/192.0);
    }
   return h_vector;
}

fn G(grids: &Grids, df_star: Vec<f64>)-> Vec<f64>{
    let vector_lenght = grids.current_grid.f.len();
    let mut g_vector = Vec::new();

    for i in 0..vector_lenght{
        g_vector.push(5.0 * (grids.previous_grid.df[i] - grids.next_grid.df[i]) / 192.0);
    }
   return g_vector;
}

// pub fn G(grids: Grids)-> Vec<f64>{
//     let vector_lenght = grids.current_grid.f.len();
//     let mut g_vector = Vec::new();

//     for i in 0..vector_lenght{
//         g_vector.push(5.0*(f_derivative(grids,i) - f_derivative(grids,i))/192.0);
//     }
//    return g_vector;
// }

// fn f_derivative(grids: Grids,index: i32)->Vec<f64>{


// }
pub fn calculate_GPlus(grids: Grids){}
// cip_internal::H(grids)
//         + cip_internal::G(grids)
//         - (cip_internal::dF(grids) - FMinus);