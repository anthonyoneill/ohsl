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

#[test]
fn test_complex_square_root() {
    let z = Cmplx::new( 1.0, 1.0 );
    let sqrt = z.sqrt();
    assert!( (sqrt.real - 1.098684113467810).abs() < 1.0e-15 );
    assert!( (sqrt.imag - 0.455089860562227).abs() < 1.0e-15 );
}

#[test]
fn test_complex_pow() {
    let z = Cmplx::new( 1.0, 1.0 );
    let w = Cmplx::new( 2.0, 1.0 );
    let pow = z.pow( &w );
    assert!( (pow.real + 0.309743504928494).abs() < 1.0e-15 );
    assert!( (pow.imag - 0.857658012588736).abs() < 1.0e-15 );
}

#[test]
fn test_complex_powf() {
    let z = Cmplx::new( 2.0, 1.0 );
    let pow = z.powf( 2.0 );
    assert!( (pow.real - 3.0).abs() < 1.0e-15 );
    assert!( (pow.imag - 4.0).abs() < 1.0e-15 );
    let pow = z.powf( 2.5 );
    assert!( (pow.real - 2.9917970717860150).abs() < 1.0e-15 );
    assert!( (pow.imag - 6.8520690100689566).abs() < 1.0e-15 );
}

#[test]
fn test_complex_exponential() {
    let z = Cmplx::new( 1.0, 1.0 );
    let exp = z.exp();
    assert!( (exp.real - 1.468693939915885).abs() < 1.0e-15 );
    assert!( (exp.imag - 2.287355287178842).abs() < 1.0e-15 );
}

#[test]
fn test_complex_logarithm() {
    let z = Cmplx::new( 1.0, 1.0 );
    let ln = z.ln();
    assert!( (ln.real - 0.346573590279973).abs() < 1.0e-15 );
    assert!( (ln.imag - 0.785398163397448).abs() < 1.0e-15 );
}

#[test]
fn test_complex_base_b_logarithm() {
    let z = Cmplx::new( 1.0, 1.0 );
    let log = z.log( Cmplx::new( 2.0, 0.0 ) );
    assert!( (log.real - 0.5).abs() < 1.0e-15 );
    assert!( (log.imag - 1.133090035456798).abs() < 1.0e-15 );
    let log = z.log( Cmplx::new( 2.0, 1.0 ) );
    assert!( (log.real - 0.745520263590820).abs() < 1.0e-15 );
    assert!( (log.imag - 0.546450996741907).abs() < 1.0e-15 );
}