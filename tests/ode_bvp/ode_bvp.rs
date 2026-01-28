use ohsl::ode_bvp::ODEBVP;
use ohsl::vector::Vec64;
use ohsl::constant::PI;
use ohsl::matrix::Mat64;

#[test]
fn constructor() {
    let nodes = Vec64::linspace( 0.0, 1.0, 5 );
    let _ode_bvp = ODEBVP::<f64>::new( nodes, 2 );
}

#[test]
fn parameters() {
    let nodes = Vec64::linspace( 0.0, 1.0, 5 );
    let mut ode_bvp = ODEBVP::<f64>::new( nodes, 2 );
    ode_bvp.tolerance( 1.0e-6 );
    ode_bvp.delta( 1.0e-7 );
    ode_bvp.iterations( 25 );  
    assert!( ode_bvp.tol == 1.0e-6 );
    assert!( ode_bvp.delta == 1.0e-7 );
    assert!( ode_bvp.max_iter == 25 );     
}

fn airyfunc( y: &Vec64, x: &f64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = y[1];
    f[1] = x * y[0];
    /*
        (y)' = y',
        (y')' = x * y,
        Airy equation y'' - x * y = 0
    */
    f
}

fn plate_bc( y: &Vec64 ) -> Vec64 {
    let mut bc = Vec64::new( 1, 0.0 );
    bc[0] = y[0];       // y(0) = 0
    bc
}

fn free_bc( y: &Vec64 ) -> Vec64 {
    let mut bc = Vec64::new( 1, 0.0 );
    bc[0] = y[0] - 1.0; // y(inf) = 1
    bc
}

#[test]
fn airy_eqn() {
    let inf = 10.0;
    let n_nodes = 20_000;
    let nodes = Vec64::linspace( 0.0, inf, n_nodes );
    let mut ode_bvp = ODEBVP::<f64>::new( nodes.clone(), 2 );
    ode_bvp.iterations( 20 );

    // Set the initial guess for the solution
    for i in 0..n_nodes {
        let x = ode_bvp.solution.coord( i );
        ode_bvp.solution[i][0] = x;     // y = x
        ode_bvp.solution[i][1] = 1.;    // y' = 1
        //TODO could use a better initial guess to match inf bc?
    }

    // Solve the BVP
    let result = ode_bvp.solve( &airyfunc, &plate_bc, &free_bc, None );
    assert!( result.is_ok() && !result.is_err() );

    // Check the solution
    let x = 8.0;
    let y = ode_bvp.solution.get_interpolated_vars( x )[0];
    let y_exact = 0.0026327428828209; 
    // ( sqrt(3) * Ai(x) - Bi(x) ) / ( sqrt(3) * Ai(x_inf) - Bi(x_inf) )
    //println!( "Airy solution at x = {}: y = {:.6e}", x, y );
    assert!( ( y - y_exact ).abs() < 1.0e-8 );
}

fn equation( y: &Vec64, _x: &f64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = y[1];
    f[1] = y[0].sin();
    /*
        (y)' = y',
        (y')' = sin( y ),
        Pendulum equation y'' - sin( y ) = 0
    */
    f
}

fn zero_bc( y: &Vec64 ) -> Vec64 {
    let mut bc = Vec64::new( 1, 0.0 );
    bc[0] = y[0];       // y(0) = 0
    bc
}

fn one_bc( y: &Vec64 ) -> Vec64 {
    let mut bc = Vec64::new( 1, 0.0 );
    bc[0] = y[0] - PI; // y(1) = pi
    bc
}

fn jacobian( y: &Vec64, _x: &f64 ) -> Mat64 {
    let mut jac = Mat64::new( 2, 2, 0.0 );
    jac[ (0, 0) ] = 0.0;            // dR1/dy1 = dR1/dy
    jac[ (0, 1) ] = 1.0;            // dR1/dy2 = dR1/dy'
    jac[ (1, 0) ] = y[0].cos();     // dR2/dy1 = dR2/dy
    jac[ (1, 1) ] = 0.0;            // dR2/dy2 = dR2/dy'
    jac
}

#[test]
fn test_eqn() {
    let n_nodes = 200;
    let nodes = Vec64::linspace( 0.0, 1.0, n_nodes );
    let mut ode_bvp = ODEBVP::<f64>::new( nodes.clone(), 2 );
    ode_bvp.iterations( 20 );

    // Set the initial guess for the solution
    for i in 0..n_nodes {
        let x = ode_bvp.solution.coord( i );
        ode_bvp.solution[i][0] = x;     // y = x
        ode_bvp.solution[i][1] = 1.;    // y' = 1
        //TODO could use a better initial guess to match inf bc?
    }

    // Solve the BVP
    let result = ode_bvp.solve( &equation, &zero_bc, &one_bc, Some(&jacobian) );
    assert!( result.is_ok() && !result.is_err() );

    // Check the solution
    let x = 1.0;
    let y = ode_bvp.solution.get_interpolated_vars( x )[0];
    let y_exact = PI;
    //println!( "Pendulum solution at x = {}: y = {:.6e}", x, y );
    assert!( ( y - y_exact ).abs() < 1.0e-8 );
}