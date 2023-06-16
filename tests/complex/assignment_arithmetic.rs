use ohsl::complex::Cmplx;

#[test]
fn add_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans += z1;
    assert_eq!( zans.real, 3.0 );
    assert_eq!( zans.imag, 3.0 );
}

#[test]
fn subtract_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans -= z1;
    assert_eq!( zans.real, 1.0 );
    assert_eq!( zans.imag, -1.0 );
}

#[test]
fn multiply_assign() {
    let z1 = Cmplx::new( 1.0, 2.0 );
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans *= z1;
    assert_eq!( zans.real, 0.0 );
    assert_eq!( zans.imag, 5.0 );
}

#[test]
fn divide_assign() {
    let mut zans = Cmplx::new( 1.0, 2.0 );
    let z2 = Cmplx::new( 2.0, 1.0 );
    zans /= z2;
    assert_eq!( zans.real, 0.8 );
    assert_eq!( zans.imag, 0.6 );
}

#[test]
fn add_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans += r;
    assert_eq!( zans.real, 6.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn subtract_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans -= r;
    assert_eq!( zans.real, -2.0 );
    assert_eq!( zans.imag, 1.0 );
}

#[test]
fn multiply_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans *= r;
    assert_eq!( zans.real, 8.0 );
    assert_eq!( zans.imag, 4.0 );
}

#[test]
fn divide_assign_real() {
    let r: f64 = 4.0;
    let mut zans = Cmplx::new( 2.0, 1.0 );
    zans /= r;
    assert_eq!( zans.real, 0.5 );
    assert_eq!( zans.imag, 0.25 );
}