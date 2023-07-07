use ohsl::complex::{Complex, Cmplx};

#[test]
fn construction() {
    let z = Complex::<i32>::new( 1, 1 );
    assert_eq!( z.real, 1 );
    assert_eq!( z.imag, 1 );
    let z = Cmplx::new( 1.0, 1.0 );
    assert_eq!( z.real, 1.0 );
    assert_eq!( z.imag, 1.0 );
}

#[test]
fn conjugate() {
    let z = Complex::new( 1.0, 1.0 );
    let zbar = z.conj();
    assert_eq!( zbar.real, 1.0 );
    assert_eq!( zbar.imag, -1.0 );
}

mod basic_arithmetic;
mod assignment_arithmetic;
mod elementary_functions;
mod trig_functions;

#[test]
fn test_complex_ordering() {
    let z = Complex::new( 0.0, 1.0 );
    let w = Complex::new( 0.0, 1.0 );
    assert!( z == w );
    let u = Complex::new( 1.0, -1.0 );
    assert!( u > z );
    let p = Complex::new( 0.0, -1.0 );
    assert!( p < z );
    assert!( u > p );
}