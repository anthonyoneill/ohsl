use ohsl::vector::Vec64;
use ohsl::matrix::Mat64;
use ohsl::newton::Newton;
use ohsl::complex::Cmplx;

#[test]
fn constructor() {
    let newton = Newton::<f64>::new( 0.0 );
    let parameters = newton.parameters();
    assert_eq!( parameters.0, 1.0e-8 );
    assert_eq!( parameters.1, 1.0e-8 );
    assert_eq!( parameters.2, 20 );
    assert_eq!( parameters.3, 0.0 );
}

#[test]
fn parameters() {
    let mut newton = Newton::<f64>::new( 0.0 );
    newton.tolerance( 1.0e-6 );
    newton.delta( 1.0e-7 );
    newton.iterations( 25 );
    newton.guess( 1.0 );
    let parameters = newton.parameters();
    assert_eq!( parameters.0, 1.0e-6 );
    assert_eq!( parameters.1, 1.0e-7 );
    assert_eq!( parameters.2, 25 );
    assert_eq!( parameters.3, 1.0 );
}

fn x_squared_minus_four(x: f64) -> f64 {
    x * x - 4.0
}

#[test]
fn solve_f64() {
    let mut newton = Newton::<f64>::new( 1.0 );
    let solution = newton.solve( &x_squared_minus_four );
    assert_eq!( solution.unwrap(), 2.0 );
    newton.guess( -1.0 );
    let solution = newton.solve( &x_squared_minus_four ).unwrap();
    assert_eq!( solution, -2.0 );
}

fn vecfunc( x: Vec64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = f64::powf( x[0], 3.0 ) + x[1] -1.0;
    f[1] = f64::powf( x[1], 3.0 ) - x[0] + 1.0;
    /*
        x^3 + y - 1 = 0,
        y^3 - x + 1 = 0,
        (x,y) = (1,0) is the only (real) solution
    */
    f
}

#[test]
fn jacobian() {
    let point = Vec64::create( vec![ 1.0, 1.0 ] );
    let jacobian = Mat64::jacobian( point, &vecfunc, 1.0e-8 );
    assert_eq!( jacobian.rows(), 2 );
    assert_eq!( jacobian.cols(), 2 );
    assert!( ( jacobian[(0,0)] - 3.0 ).abs() < 1.0e-6 );
    assert!( ( jacobian[(0,1)] - 1.0 ).abs() < 1.0e-6 );
    assert!( ( jacobian[(1,0)] + 1.0 ).abs() < 1.0e-6 );
    assert!( ( jacobian[(1,1)] - 3.0 ).abs() < 1.0e-6 );
}

#[test]
fn solve_vec64() {
    let guess = Vec64::create( vec![ 0.5, 0.25 ] );
    let newton = Newton::<Vec64>::new( guess );
    let solution = newton.solve( &vecfunc ).unwrap();
    assert!( ( solution[0] - 1.0) < 1.0e-8 );
    assert!( ( solution[1] - 0.0) < 1.0e-8 );
}

fn jacfunc( x: Vec64 ) -> Mat64 {
    let mut j = Mat64::new( 2, 2, 0.0 );
    j[(0,0)] = 3.0 * f64::powf( x[0], 2.0 );
    j[(0,1)] = 1.0;
    j[(1,0)] = -1.0;
    j[(1,1)] = 3.0 * f64::powf( x[1], 2.0 );
    /*
        J = [ 3x^2,  1   
                -1 , 3y^2 ]
    */
    j
}

#[test]
fn exact_jacobian_solve() {
    let guess = Vec64::create( vec![ 0.5, 0.25 ] );
    let newton = Newton::<Vec64>::new( guess );
    let solution = newton.solve_jacobian( &vecfunc, &jacfunc ).unwrap();
    assert!( ( solution[0] - 1.0) < 1.0e-8 );
    assert!( ( solution[1] - 0.0) < 1.0e-8 );
}

fn x_cubed_minus_two(x: Cmplx) -> Cmplx {
    x * x * x - 2.0
}

#[test]
fn solve_cmplx() {
    let mut newton = Newton::<Cmplx>::new( Cmplx::new( 1.0, 0.0 ) );
    let solution = newton.solve( &x_cubed_minus_two );
    assert_eq!( solution.unwrap(), Cmplx::new( 1.2599210498948732, 0.0 ) );
    newton.guess( Cmplx::new( -1.0, 1.0 ) );
    let solution = newton.solve( &x_cubed_minus_two );
    assert_eq!( solution.unwrap(), Cmplx::new( -0.6299605249474366, 1.0911236359717214 ) );
    newton.guess( Cmplx::new( -1.0, -1.0 ) );
    let solution = newton.solve( &x_cubed_minus_two );
    assert_eq!( solution.unwrap(), Cmplx::new( -0.6299605249474366, -1.0911236359717214 ) );
}