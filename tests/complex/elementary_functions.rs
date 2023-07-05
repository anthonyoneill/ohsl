use ohsl::complex::{Complex, Cmplx};
use ohsl::traits::{Zero, One};

#[test]
fn test_complex_abs_sqr() {
    let z = Complex::new( 1.0, 1.0 );
    let r = z.abs_sqr();
    assert_eq!( r, 2.0 );
}

#[test]
fn test_complex_abs() {
    let z = Complex::new( 3.0, 4.0 );
    let r = z.abs();
    assert_eq!( r, 5.0 );
}

#[test]
fn test_complex_arg() {
    let z = Cmplx::new( 3.0, 4.0 );
    let theta = z.arg();
    assert_eq!( theta, f64::atan( 4.0 / 3.0 ) );
}

#[test]
fn test_complex_zero() {
    let z = Complex::<i32>::zero();
    assert_eq!( z.real, 0 );
    assert_eq!( z.imag, 0 );
}

#[test]
fn test_complex_one() {
    let z = Complex::<i32>::one();
    assert_eq!( z.real, 1 );
    assert_eq!( z.imag, 0 );
}