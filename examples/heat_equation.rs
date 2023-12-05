// cargo run --example heat_equation
use std::fs;
use ohsl::{Vec64, Mesh2D, Tridiagonal};

fn main() {
    println!("----- Heat equation -----");

    const T_MAX: f64 = 1.0;         // Maximum simulation time
    const N: usize = 200;           // Number of time steps
    const J: usize = 50;            // Number of space steps
    let dt = T_MAX / N as f64;      // Temporal step size
    let dx = 1.0 / J as f64;        // Spatial step size
    let alpha = 1.0;                // Thermal diffusion coefficient

    println!( " * Solving the heat equation u_t = \u{03B1} u_xx in the domain" );
    println!( " * t=[0,{}] and x=[0,1] subject to the boundary conditions", T_MAX );
    println!( " * u(0,t) = u(1,t) = 0 and the initial condition u(x,0) = 100 * sin( pi * x ).");
    println!( " * The time step size is dt = {} and the spacial step size is dx = {}.", dt, dx );
    println!( " * The number of time and space steps are N = {} and J = {} respectively.", N, J);
    println!( " * The thermal diffusion coefficient is \u{03B1} = {}.", alpha );

    let mu = ( alpha * dt ) / ( dx * dx ); 

    // Mesh for storing the solution
    let tnodes = Vec64::linspace( 0.0, T_MAX, N + 1 );
    let xnodes = Vec64::linspace( 0.0, 1.0, J + 1 );
    let mut solution = Mesh2D::<f64>::new( xnodes, tnodes, 1 );

    // Matrices for Crank-Nicolson method including BCs u(0,t) = u(1,t) = 0
    let mut implicit = Tridiagonal::with_elements( -mu, 2. + 2. * mu, -mu, J + 1 );
    implicit[( 0, 0 )]     = 1.0;
    implicit[( 0, 1 )]     = 0.0;
    implicit[( J, J - 1 )] = 0.0;
    implicit[( J, J )]     = 1.0;

    let mut explicit = Tridiagonal::with_elements( mu, 2. - 2. * mu, mu, J + 1 );
    explicit[( 0, 0 )]     =   1.0;
    explicit[( 0, 1 )]     =   0.0;
    explicit[( J, J - 1 )] =   0.0;
    explicit[( J, J )]     = - 1.0;

    // Initial condition u(x,0) = 100 * sin( pi * x )
    let mut current = Vec64::zeros( J + 1 );
    current[ 0 ] = 0.0;
    solution[( 0, 0 )][ 0 ] = current[ 0 ];
    for j in 1..J {
        current[ j ] = 100.0 * ( std::f64::consts::PI * ( j as f64 ) * dx ).sin();
        solution[( j, 0 )][ 0 ] = current[ j ];
    }
    current[ J ] = 0.0;
    solution[( J, 0 )][ 0 ] = current[ J ];
    let mut next: Vec64;

    println!( " * Solving the heat equation using the Crank-Nicolson method..." );
    // Time loop
    for i in 1..N + 1 {
        current = &explicit * &current;
        next = implicit.solve( &current );
        current = next;
        for j in 0..J + 1 {
            solution[( j, i )][ 0 ] = current[ j ];
        }
    }
    println!( " * Finished solving the heat equation." );
    fs::create_dir_all( "DATA" ).expect( "Could not create DATA output directory" );
    let mut string = String::from( "DATA/Heat_equation" );
    string.push_str( format!( "_J_{}_N_{}.dat", J, N ).as_str() );
    solution.output( &string, 6 );

    print!( " * For an animation of the solution run: " );
    println!( "python3 Plotting/Heat_plot.py {} {}", J, N );
    println!( " * and play the resulting file DATA/Heat_equation.avi" );
    
    println!( "--- FINISHED ---" );
}