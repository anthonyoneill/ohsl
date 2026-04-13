use ohsl::ode_ivp::ODEIVP;
use ohsl::vector::Vec64;

fn lorenz( y: &Vec64, _x: &f64 ) -> Vec64 {
    let a = 10.0;
    let b = 7.0;
    let c = 8.0 / 3.0;
    let mut f = Vec64::new( 3, 0.0 );
    f[0] = a * ( y[1] - y[0] );              
    f[1] = y[0] * ( b - y[2] ) - y[1];             
    f[2] = y[0] * y[1] - c * y[2];
    /*
        Lorenz system:
        x' = a * (y - x)
        y' = x * (b - z) - y
        z' = x * y - c * z
    */
    f
}

fn main() {
    println!("----- Lorenz system -----");

    println!( " * Solving the Lorenz system" );
    println!( "\t x' = a * (y - x), ");
    println!( "\t y' = x * (b - z) - y, ");
    println!( "\t z' = x * y - c * z, ");
    println!( " * in the time domain t=[0,200] with initial conditions x(0)=y(0)=z(0)=1." );
    println!( " * The parameters are chosen to be a=10, b=7, c=8/3." );
    println!( " * There is a fixed point solution at (x,y,z) = (4,4,6)." );

    let max_steps = 5_000_000;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 200.0, max_steps );
    let initial_conditions = Vec64::new( 3, 1.0 );
    let result = ode_ivp.rkf45( &lorenz, initial_conditions );
    match result {
        Ok( () ) => {
            println!( " * Successfully solved the Lorenz system with RKF45 method." );
            println!( " * The solution has {} nodes.", ode_ivp.solution.nodes().size() );
        },
        Err( msg ) => {
            println!( " * Error solving the Lorenz system: {}", msg );
        }
    }
    let sol_n = ode_ivp.solution.get_interpolated_vars( 200.0 );
    println!( " * At t = 200: x = {:.6}, y = {:.6}, z = {:.6}", sol_n[0], sol_n[1], sol_n[2] );
    let exact = Vec64::create( vec![ 4.0, 4.0, 6.0 ] );
    let error = ( sol_n - exact ).norm_inf();
    println!( " * Error at t = 200 compared to reference solution: {:.2e}", error );
    println!( " * Finished solving the Lorenz system." );
}