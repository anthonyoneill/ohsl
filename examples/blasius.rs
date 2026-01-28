// cargo run --example blasius
use std::fs;
use ohsl::ode_bvp::ODEBVP;
use ohsl::vector::{Vec64, Vector};

fn blasius_eqn( y: &Vec64, _x: &f64 ) -> Vec64 {
    let mut f = Vec64::new( 3, 0.0 );
    f[0] = y[1];
    f[1] = y[2];
    f[2] = -y[0] * y[2];
    /*
        (y)' = y',
        (y')' = y'',
        (y'')' = - y * y'' ,
        Blasius equation y''' + y * y'' = 0
    */
    f
}

fn plate_bc( y: &Vec64 ) -> Vector<Option<f64>> {
    let mut bc = Vector::<Option<f64>>::new( 3, None );
    bc[0] = Some( y[0] );       // y(0) = 0
    bc[1] = Some( y[1] );       // y'(0) = 0
    bc[2] = None;               // y''(0) free
    bc
}

fn free_bc( y: &Vec64 ) -> Vector<Option<f64>> {
    let mut bc = Vector::<Option<f64>>::new( 3, None );
    bc[0] = None;               // y(0) free
    bc[1] = Some( y[1] - 1.0 ); // y'(inf) = 1
    bc[2] = None;               // y''(inf) free
    bc
}

fn main() {
    println!("----- Blasius equation -----");

    const ETA_INF: f64 = 20.0;
    const N: usize = 512;

    println!( " * Solving the Blasius equation f''' + f * f'' = 0 in the domain" );
    println!( " * \u{03B7}=[0,{}] subject to the boundary conditions", ETA_INF );
    println!( " * f(0) = 0, f'(0) = 0, f' -> 1 as \u{03B7} -> \u{03B7}_\u{221E} = {}.", ETA_INF );
    println!( " * The number of nodes in the discretisation is N = {}.", N );

    let nodes = Vec64::powspace( 0.0, ETA_INF, N, 1.2 ); // More nodes near eta=0
    let mut ode_bvp = ODEBVP::<f64>::new( nodes.clone(), 3 );
    ode_bvp.iterations( 20 );

    // Set the initial guess for the solution
    for i in 0..N {
        let eta = ode_bvp.solution.coord( i );
        let exp_minus_eta = (-eta).exp();
        ode_bvp.solution[i][0] = eta * ( 1.0 - exp_minus_eta ); 
        ode_bvp.solution[i][1] = 1.0 - exp_minus_eta + eta * exp_minus_eta; 
        ode_bvp.solution[i][2] = 2.0 * exp_minus_eta - eta * exp_minus_eta;
    }

    println!( " * Solving the Blasius equation using the finite difference method..." );
    let result = ode_bvp.solve( &blasius_eqn, &plate_bc, &free_bc, None );  
    assert!( result.is_ok() && !result.is_err() );

    let fdd_0 = ode_bvp.solution.get_interpolated_vars( 0.0 )[2];
    println!( " * The computed value of f''(0) is {:.8}.", fdd_0 ); 

    let c= 1.65519036023f64;
    let reference = 1. / c.powf( 1.5 );
    println!( " * The reference value of f''(0) is {:.8}.", reference );
    println!( " * The relative error is {:.8}%.", ( fdd_0 - reference ).abs() / reference * 100.0 );  

    println!( " * Finished solving the Blasius equation." );
    fs::create_dir_all( "DATA" ).expect( "Could not create DATA output directory" );
    let mut string = String::from( "DATA/Blasius_equation" );
    string.push_str( format!( "_N_{}.dat", N ).as_str() );
    ode_bvp.solution.output( &string, 6 ); 

    print!( " * For a plot of the solution run: " );
    println!( "python3 Plotting/Blasius_plot.py {}", N );
    println!( " * and view the resulting file DATA/Blasius_equation_N_{}.png", N );
    
    println!( "--- FINISHED ---" );
}