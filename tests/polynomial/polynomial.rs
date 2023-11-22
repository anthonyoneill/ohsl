use ohsl::{Polynomial, Cmplx};

#[test]
fn unspecified_polynomial() {
    let p = Polynomial::<f64>::empty();
    assert_eq!( p.size(), 0 );
}

#[test]
fn specified_polynomial() {
    let mut p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] );
    assert_eq!( p.size(), 3 );
    assert_eq!( p[0], 1.0 );
    assert_eq!( p[1], 2.0 );
    assert_eq!( p[2], 3.0 );
    p[2] = 4.0;
    assert_eq!( p[2], 4.0 );
    assert_eq!( p.degree(), 2 );
    assert_eq!( p.coeffs(), &vec![ 1.0, 2.0, 4.0 ] );
    *p.coeffs() = vec![ 1.0, 1.0, 1.0 ];
    assert_eq!( p.coeffs(), &vec![ 1.0, 1.0, 1.0 ] );
}

#[test]
fn evaluate_polynomial() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    assert_eq!( p.evaluate( 0.0 ), 1.0 );
    assert_eq!( p.evaluate( 1.0 ), 6.0 );
    assert_eq!( p.evaluate( 2.0 ), 17.0 );
}

#[test]
fn quadratic_polynomial() {
    let mut p = Polynomial::quadratic( 1.0, 2.0, 4.0 ); // x^2 + 2x + 4
    assert_eq!( p.size(), 3 );
    assert_eq!( p.degree(), 2 );
    assert_eq!( p[0], 4.0 );
    assert_eq!( p[1], 2.0 );
    assert_eq!( p[2], 1.0 );
    assert_eq!( p.coeffs(), &vec![ 4.0, 2.0, 1.0 ] );
    assert_eq!( p.evaluate( 2.0 ), 12.0 );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 2 );
    let r0 = Cmplx::new( -1.0, -(3.0_f64).sqrt() );
    let r1 = Cmplx::new( -1.0,  (3.0_f64).sqrt() );
    let eps = 1.0e-15; //TODO use f64::EPSILON when refinement is implemented
    assert!( ( roots[0].real - r0.real ).abs() < eps );
    assert!( ( roots[0].imag - r0.imag ).abs() < eps );
    assert!( ( roots[1].real - r1.real ).abs() < eps );
    assert!( ( roots[1].imag - r1.imag ).abs() < eps );
}

#[test]
fn cubic_polynomial() {
    let mut p = Polynomial::cubic( 1.0, 2.0, 3.0, 4.0 ); // x^3 + 2x^2 + 3x + 4
    assert_eq!( p.size(), 4 );
    assert_eq!( p.degree(), 3 );
    assert_eq!( p[0], 4.0 );
    assert_eq!( p[1], 3.0 );
    assert_eq!( p[2], 2.0 );
    assert_eq!( p[3], 1.0 );
    assert_eq!( p.coeffs(), &vec![ 4.0, 3.0, 2.0, 1.0 ] );
    assert_eq!( p.evaluate( 1.0 ), 10.0 );
    assert_eq!( p.evaluate( 2.0 ), 26.0 );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 3 );
    let r0 = Cmplx::new( -1.65062919143939, 0.0 );
    let r1 = Cmplx::new( -0.17468540428031, -1.54686888723140 );
    let r2 = Cmplx::new( -0.17468540428031,  1.54686888723140 );
    let eps = 1.0e-14; //TODO use f64::EPSILON when refinement is implemented
    assert!( ( roots[0].real - r0.real ).abs() < eps );
    assert!( ( roots[0].imag - r0.imag ).abs() < eps );
    assert!( ( roots[1].real - r1.real ).abs() < eps );
    assert!( ( roots[1].imag - r1.imag ).abs() < eps );
    assert!( ( roots[2].real - r2.real ).abs() < eps );
    assert!( ( roots[2].imag - r2.imag ).abs() < eps );
}

#[test]
fn cubic_repeated_roots() {
    let p = Polynomial::cubic( 1.0, 3.0, 3.0, 1.0 ); // x^3 + 3x^2 + 3x + 1
    let roots = p.roots( true );
    assert_eq!( roots.size(), 3 );
    let r0 = Cmplx::new( -1.0, 0.0 );
    let r1 = Cmplx::new( -1.0, 0.0 );
    let r2 = Cmplx::new( -1.0, 0.0 );
    let eps = 1.0e-14; //TODO use f64::EPSILON when refinement is implemented
    assert_eq!( roots[0], r0 );
    assert!( ( roots[0].real - r0.real ).abs() < eps );
    assert!( ( roots[0].imag - r0.imag ).abs() < eps );
    assert!( ( roots[1].real - r1.real ).abs() < eps );
    assert!( ( roots[1].imag - r1.imag ).abs() < eps );
    assert!( ( roots[2].real - r2.real ).abs() < eps );
    assert!( ( roots[2].imag - r2.imag ).abs() < eps );
}

//TODO need more testing of cubic_solve