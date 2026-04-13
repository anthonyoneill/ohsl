use ohsl::ode_ivp::ODEIVP;
use ohsl::vector::Vec64;

#[test]
fn constructor() {
    let ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, 5 );
    assert_eq!( ode_ivp.x0, 0.0 );
    assert_eq!( ode_ivp.xn, 1.0 );
    assert_eq!( ode_ivp.h, 0.2 );
    assert_eq!( ode_ivp.max_steps, 5 );
    assert_eq!( ode_ivp.store_every, 1 );
}

fn first_order_func( y: &Vec64, x: &f64 ) -> Vec64 {
    let mut f = Vec64::new( 1, 0.0 );
    f[0] = 3.0 * x.exp() - 2.0 * y[0];
    /*
        (y)' = 3e^x - 2y,
    */
    f
}

#[test]
fn rk4_first_order() {
    let n = 20;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, n );
    let initial_conditions = Vec64::new( 1, 3.0 ); // y(0) = 3
    let result = ode_ivp.rk4( &first_order_func, initial_conditions );
    assert!( result.is_ok() );
    //let y_final = ode_ivp.solution.get_nodes_vars( n );
    let y_final = ode_ivp.solution.get_interpolated_vars( 1.0 );
    let y_exact = 2.0 * (-2.0f64).exp() + 1.0f64.exp();
    let error = ( y_final[0] - y_exact ).abs();
    /*println!( 
        "RK4 solution at x = 1: y = {:.6e}, exact = {:.6e}, error = {:.2e}", 
        y_final[0], y_exact, error 
    );*/
    assert!( error < 1.0e-4 );
}

#[test]
fn store_every() {
    let n = 21;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, n );
    ode_ivp.store_every( 5 ); // Store every 5 steps
    let initial_conditions = Vec64::new( 1, 3.0 ); // y(0) = 3
    let result = ode_ivp.rk4( &first_order_func, initial_conditions );
    assert!( result.is_ok() );
    assert_eq!( ode_ivp.solution.nodes().size(), 6 ); 
    let y_final = ode_ivp.solution.get_interpolated_vars( 1.0 );
    let y_exact = 2.0 * (-2.0f64).exp() + 1.0f64.exp();
    let error = ( y_final[0] - y_exact ).abs();
    /*println!( 
        "RK4 solution at x = 1: y = {:.6e}, exact = {:.6e}, error = {:.2e}", 
        y_final[0], y_exact, error 
    );*/
    assert!( error < 1.0e-4 );
}

fn second_order_func( y: &Vec64, _x: &f64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = y[1];              
    f[1] = -y[0];             
    /*
        (y)' = y',
        (y')' = -y,
        Simple harmonic oscillator
    */
    f
}

#[test]
fn rk4_second_order() {
    let n = 100;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 2.0 * std::f64::consts::PI, n );
    let mut initial_conditions = Vec64::new( 2, 0.0 ); // y(0) = 0, y'(0) = 1
    initial_conditions[1] = 1.0; // y'(0) = 1
    let result = ode_ivp.rk4( &second_order_func, initial_conditions );
    assert!( result.is_ok() );
    let y_final = ode_ivp.solution.get_interpolated_vars( 2.0 * std::f64::consts::PI );
    let mut y_exact = Vec64::new( 2, 0.0 );
    y_exact[0] = 0.0; // sin(2pi) = 0
    y_exact[1] = 1.0; // cos(2pi) = 1
    let error = ( y_final.clone() - y_exact.clone() ).norm_inf();
    /*println!( "RK4 solution at x = 2pi: y = [{:.6e}, {:.6e}]", y_final[0], y_final[1] );
    println!( "Exact solution at x = 2pi: y = [{:.6e}, {:.6e}]", y_exact[0], y_exact[1] );
    println!( "Error at x = 2pi: {:.2e}", error );*/
    assert!( error < 1.0e-4 );
} 

#[test]
fn rk4_incorrect_initial_conditions() {
    let n = 20;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, n );
    let initial_conditions = Vec64::new( 2, 3.0 ); // Incorrect size
    let result = ode_ivp.rk4( &first_order_func, initial_conditions );
    assert!( result.is_err() );
}

#[test]
fn rkf45_first_order() {
    let max_steps = 100;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, max_steps );
    let initial_conditions = Vec64::new( 1, 3.0 ); // y(0) = 3
    let result = ode_ivp.rkf45( &first_order_func, initial_conditions );
    match result {
        Ok( () ) => println!( "RKF45 completed successfully." ),
        Err( ref e ) => println!( "RKF45 failed: {}", e ),
    }
    assert!( result.is_ok() );
    //let y_final = ode_ivp.solution.get_nodes_vars( n );
    let y_final = ode_ivp.solution.get_interpolated_vars( 1.0 );
    let y_exact = 2.0 * (-2.0f64).exp() + 1.0f64.exp();
    let error = ( y_final[0] - y_exact ).abs();
    /*println!( 
        "RKF45 solution at x = 1: y = {:.6e}, exact = {:.6e}, error = {:.2e}", 
        y_final[0], y_exact, error 
    );*/
    assert!( error < ode_ivp.tol );

    let num_nodes = ode_ivp.solution.nodes().size();
    println!( "RKF45 number of nodes: {}", num_nodes );
    assert!( num_nodes < max_steps ); // Should not exceed the maximum number of steps
}

#[test]
fn rkf45_incorrect_initial_conditions() {
    let max_steps = 100;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 1.0, max_steps );
    let initial_conditions = Vec64::new( 2, 3.0 ); // Incorrect size
    let result = ode_ivp.rkf45( &first_order_func, initial_conditions );
    assert!( result.is_err() );
}

#[test]
fn rkf45_second_order() {
    let max_steps = 100;
    let mut ode_ivp = ODEIVP::<f64>::new( 0.0, 2.0 * std::f64::consts::PI, max_steps );
    ode_ivp.tolerance( 1.0e-6 );
    let mut initial_conditions = Vec64::new( 2, 0.0 ); // y(0) = 0, y'(0) = 1
    initial_conditions[1] = 1.0; // y'(0) = 1
    let result = ode_ivp.rkf45( &second_order_func, initial_conditions );
    match result {
        Ok( () ) => println!( "RKF45 completed successfully." ),
        Err( ref e ) => println!( "RKF45 failed: {}", e ),
    }
    assert!( result.is_ok() );
    let y_final = ode_ivp.solution.get_interpolated_vars( 2.0 * std::f64::consts::PI );
    let mut y_exact = Vec64::new( 2, 0.0 );
    y_exact[0] = 0.0; // sin(2pi) = 0
    y_exact[1] = 1.0; // cos(2pi) = 1
    let error = ( y_final.clone() - y_exact.clone() ).norm_inf();
    /*println!( "RKF45 solution at x = 2pi: y = [{:.6e}, {:.6e}]", y_final[0], y_final[1] );
    println!( "Exact solution at x = 2pi: y = [{:.6e}, {:.6e}]", y_exact[   0], y_exact[1] );
    println!( "Error at x = 2pi: {:.2e}", error );*/
    assert!( error < ode_ivp.tol );
}