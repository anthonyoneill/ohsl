use ohsl::complex::{Complex, Cmplx};
use ohsl::traits::{Zero, One};

#[test]
fn test_complex_construction() {
    // Integer construction
    let z = Complex::<i32>::new( 1, 1 );
    assert_eq!( z.real, 1 );
    assert_eq!( z.imag, 1 );
    // Float construction
    let z = Cmplx::new( 1.0, 1.0 );
    assert_eq!( z.real, 1.0 );
    assert_eq!( z.imag, 1.0 );
}

#[test]
fn test_complex_conj() {
    let z = Complex::new( 1.0, 1.0 );
    let zbar = z.conj();
    assert_eq!( zbar.real, 1.0 );
    assert_eq!( zbar.imag, -1.0 );
}

#[test]
fn test_complex_unary_sub() {
    let z = Cmplx::new( 1.0, 1.0 );
    let zans = - z;
    assert_eq!( zans.real, -1.0 );
    assert_eq!( zans.imag, -1.0 );
}

#[test]
fn test_complex_binary_add() {
    let z1 = Complex::<f32>::new( 1.0, 2.0 );
    let z2 = Complex::<f32>::new( 2.0, 1.0 );
    let zans = z1 + z2;
    assert_eq!( zans.real, 3.0 );
    assert_eq!( zans.imag, 3.0 );
}

#[test]
fn test_complex_binary_sub() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 - z2;
    assert_eq!( zans.real, -1.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn test_complex_mul() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 * z2;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 5.0 );
}

#[test]
fn test_complex_div() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 / z2;
    assert_eq!( zans.real, 0.8 );
    assert_eq!( zans.imag, 0.6 );
}

#[test]
fn test_complex_plus_real() {
    let r: f64 = 1.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z + r;
    assert_eq!( zans.real, 2.0 );
    assert_eq!( zans.imag, 2.0 )
}

#[test]
fn test_complex_minus_real() {
    let r: f64 = 1.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z - r;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 2.0 )
}

#[test]
fn test_complex_times_real() {
    let r: f64 = 4.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z * r;
    assert_eq!( zans.real, 4.0 );
    assert_eq!( zans.imag, 8.0 )
}

#[test]
fn test_complex_divide_real() {
    let r: f64 = 4.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z / r;
    assert_eq!( zans.real, 0.25 );
    assert_eq!( zans.imag, 0.5 )
}

#[test]
fn test_complex_add_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans += z1;
    assert_eq!( zans.real, 3.0 );
    assert_eq!( zans.imag, 3.0 );
}

#[test]
fn test_complex_sub_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans -= z1;
    assert_eq!( zans.real, 1.0 );
    assert_eq!( zans.imag, -1.0 );
}

#[test]
fn test_complex_mul_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans *= z1;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 5.0 );
}

#[test]
fn test_complex_div_assign() {
    let mut zans = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    zans /= z2;
    assert_eq!( zans.real, 0.8 );
    assert_eq!( zans.imag, 0.6 );
}

#[test]
fn test_complex_add_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans += r;
    assert_eq!( zans.real, 6.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn test_complex_sub_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans -= r;
    assert_eq!( zans.real, -2.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn test_complex_mul_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans *= r;
    assert_eq!( zans.real, 8.0 );
    assert_eq!( zans.imag, 4.0 );
}

#[test]
fn test_complex_div_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans /= r;
    assert_eq!( zans.real, 0.5 );
    assert_eq!( zans.imag, 0.25 );
}

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
    assert_eq!( theta, libm::atan( 4.0 / 3.0 ) );
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