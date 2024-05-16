use ohsl::matrix::{Matrix, Mat64};
use ohsl::vector::{Vector, Vec64};
use ohsl::complex::Cmplx;

#[test]
fn solve_basic() {
    let mut a = Mat64::new( 2, 2, 1.0 );
    a[(0,1)] = 2.0;
    a[(1,0)] = 3.0; 
    a[(1,1)] = 4.0;
    let b = Vec64::create( vec![ 5.0, 11.0 ] ); 
    let x = a.solve_basic( &b );
    assert_eq!( x[0], 1.0 );
    assert_eq!( x[1], 2.0 );
    let mut c = Matrix::<Cmplx>::new( 2, 2, Cmplx::new( 1.0, 1.0 ) );
    c[(0,1)] = Cmplx::new( -1.0, 0.0 );
    c[(1,0)] = Cmplx::new( 1.0, -1.0 );
    let mut d = Vector::<Cmplx>::new( 2, Cmplx::new( 0.0, 1.0 ) );
    d[1] = Cmplx::new( 1.0, 0.0 );
    let y = c.solve_basic( &d );
    assert_eq!( y[0], Cmplx::new( 0.5, 0.5 ) );
    assert_eq!( y[1], Cmplx::new( 0.0, 0.0 ) );
}

#[test]
fn solve_lu() {
    let mut a = Mat64::new( 2, 2, 4.0 );
    a[(0,1)] = 3.0;
    a[(1,0)] = 6.0; 
    a[(1,1)] = 3.0;
    //let (pivots, P) = a.lu_decomp_in_place();
    let b = Vec64::create( vec![ 10.0, 12.0 ] );
    let x = a.solve_lu( &b );
    assert_eq!( x[0], 1.0 );
    assert_eq!( x[1], 2.0 );
}