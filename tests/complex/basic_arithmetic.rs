use ohsl::complex::{Complex, Cmplx};

#[test]
fn unary_subtract() {
    let z = Cmplx::new( 1.0, 1.0 );
    let zans = - z;
    assert_eq!( zans.real, -1.0 );
    assert_eq!( zans.imag, -1.0 );
}

#[test]
fn binary_add() {
    let z1 = Complex::<f32>::new( 1.0, 2.0 );
    let z2 = Complex::<f32>::new( 2.0, 1.0 );
    let zans = z1 + z2;
    assert_eq!( zans.real, 3.0 );
    assert_eq!( zans.imag, 3.0 );
}

#[test]
fn binary_subtract() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 - z2;
    assert_eq!( zans.real, -1.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn multiply() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 * z2;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 5.0 );
}

#[test]
fn divide() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    let zans = z1 / z2;
    assert_eq!( zans.real, 0.8 );
    assert_eq!( zans.imag, 0.6 );
}

#[test]
fn plus_real() {
    let r: f64 = 1.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z + r;
    assert_eq!( zans.real, 2.0 );
    assert_eq!( zans.imag, 2.0 );
}

#[test]
fn minus_real() {
    let r: f64 = 1.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z - r;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 2.0 );
}

#[test]
fn multiply_real() {
    let r: f64 = 4.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z * r;
    assert_eq!( zans.real, 4.0 );
    assert_eq!( zans.imag, 8.0 );
    let zans = r * z;
    assert_eq!( zans.real, 4.0 );
    assert_eq!( zans.imag, 8.0 );
}

#[test]
fn divide_real() {
    let r: f64 = 4.0;
    let z = Cmplx::new( 1.0, 2.0 );
    let zans = z / r;
    assert_eq!( zans.real, 0.25 );
    assert_eq!( zans.imag, 0.5 );
}