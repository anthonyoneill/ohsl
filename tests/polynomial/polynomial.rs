use std::vec;

use ohsl::{Polynomial, Cmplx, One};

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
    assert_eq!( p.degree().unwrap(), 2 );
    assert_eq!( p.coeffs(), &vec![ 1.0, 2.0, 4.0 ] );
    *p.coeffs() = vec![ 1.0, 1.0, 1.0 ];
    assert_eq!( p.coeffs(), &vec![ 1.0, 1.0, 1.0 ] );
}

#[test]
fn evaluate_polynomial() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    assert_eq!( p.eval( 0.0 ), 1.0 );
    assert_eq!( p.eval( 1.0 ), 6.0 );
    assert_eq!( p.eval( 2.0 ), 17.0 );
}

#[test]
fn quadratic_polynomial() {
    let mut p = Polynomial::quadratic( 1.0, 2.0, 4.0 ); // x^2 + 2x + 4
    assert_eq!( p.size(), 3 );
    assert_eq!( p.degree().unwrap(), 2 );
    assert_eq!( p[0], 4.0 );
    assert_eq!( p[1], 2.0 );
    assert_eq!( p[2], 1.0 );
    assert_eq!( p.coeffs(), &vec![ 4.0, 2.0, 1.0 ] );
    assert_eq!( p.eval( 2.0 ), 12.0 );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 2 );
    let r0 = Cmplx::new( -1.0, -1.732050807568877293527446 );
    let r1 = Cmplx::new( -1.0,  1.732050807568877293527446 );
    let eps = 2.0 * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
}

#[test]
fn cubic_polynomial() {
    let mut p = Polynomial::cubic( 1.0, 2.0, 3.0, 4.0 ); // x^3 + 2x^2 + 3x + 4
    assert_eq!( p.size(), 4 );
    assert_eq!( p.degree().unwrap(), 3 );
    assert_eq!( p[0], 4.0 );
    assert_eq!( p[1], 3.0 );
    assert_eq!( p[2], 2.0 );
    assert_eq!( p[3], 1.0 );
    assert_eq!( p.coeffs(), &vec![ 4.0, 3.0, 2.0, 1.0 ] );
    assert_eq!( p.eval( 1.0 ), 10.0 );
    assert_eq!( p.eval( 2.0 ), 26.0 );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 3 );
    let r0 = Cmplx::new( -1.650629191439388218880801, 0.0 );
    let r1 = Cmplx::new( -0.174685404280305890559600, -1.546868887231396277142806 );
    let r2 = Cmplx::new( -0.174685404280305890559600,  1.546868887231396277142806 );
    let eps = 2. * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
    assert!( ( roots[2] - r2 ).abs() < eps );
}

#[test]
fn cubic_repeated_roots() {
    let p = Polynomial::cubic( 1.0, 3.0, 3.0, 1.0 ); // x^3 + 3x^2 + 3x + 1
    let roots = p.roots( true );
    assert_eq!( roots.size(), 3 );
    let r0 = Cmplx::new( -1.0, 0.0 );
    let r1 = Cmplx::new( -1.0, 0.0 );
    let r2 = Cmplx::new( -1.0, 0.0 );
    let eps = 2. * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
    assert!( ( roots[2] - r2 ).abs() < eps );
}

//TODO need more testing of cubic_solve

#[test]
fn quartic_polynomial() {
    let p = Polynomial::new( vec![ 5.0, 4.0, 3.0, 2.0, 1.0 ] ); // x^4 + 2x^3 + 3x^2 + 4x + 5
    assert_eq!( p.degree().unwrap(), 4 );
    assert_eq!( p.eval( 1.0 ), 15.0 );
    assert_eq!( p.eval( -1.0 ), 3.0 );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 4 );
    let r0 = Cmplx::new( -1.287815479557647988872498, -0.8578967583284902864164198 );
    let r1 = Cmplx::new( -1.287815479557647988872498,  0.8578967583284902864164198 );
    let r2 = Cmplx::new(  0.287815479557647988872498, -1.416093080171907938724658 );
    let r3 = Cmplx::new(  0.287815479557647988872498,  1.416093080171907938724658 );
    let eps = 2. * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
    assert!( ( roots[2] - r2 ).abs() < eps );
    assert!( ( roots[3] - r3 ).abs() < eps );
}

#[test]
fn complex_quadratic() { // (1+i)x^2 + x + 1
    let p = Polynomial::quadratic( Cmplx::new( 1.0, 1.0 ), Cmplx::one(), Cmplx::one() );
    assert_eq!( p.degree().unwrap(), 2 );  
    assert_eq!( p.eval( Cmplx::one() ), Cmplx::new( 3.0, 1.0 ) );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 2 );
    let r0 = Cmplx::new( 0.0, 1.0 );
    let r1 = Cmplx::new( -0.5, -0.5 );
    let eps = 2. * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
}

#[test]
fn complex_cubic() { // (1+i)x^3 + x^2 + x + 1
    let p = Polynomial::cubic( Cmplx::new( 1.0, 1.0 ), Cmplx::one(), Cmplx::one(), Cmplx::one() );
    assert_eq!( p.degree().unwrap(), 3 );
    assert_eq!( p.eval( Cmplx::one() ), Cmplx::new( 4.0, 1.0 ) );
    let roots = p.roots( true );
    assert_eq!( roots.size(), 3 );
    let r0 = Cmplx::new( -0.7744803530832973567204905,  0.2599691354027557701265485 );
    let r1 = Cmplx::new( -0.0932651775266529118468203, -0.7875537459442802686603073 );
    let r2 = Cmplx::new(  0.3677455306099502685673108,  1.0275846105415244985337588 );
    let eps = 2. * f64::EPSILON;
    assert!( ( roots[0] - r0 ).abs() < eps );
    assert!( ( roots[1] - r1 ).abs() < eps );
    assert!( ( roots[2] - r2 ).abs() < eps );
}

#[test]
fn polydiv_divide_by_empty() {
    let u = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let v = Polynomial::empty();
    let result = u.polydiv( &v );
    assert!( result.is_err() );
}

#[test]
fn polydiv_divide_by_zero() {
    let u = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let v = Polynomial::new( vec![ 0.0, 0.0, 0.0 ] );
    let result = u.polydiv( &v );
    assert!( result.is_err() );
}

#[test]
fn polydiv_divide_by_one() {
    let u = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let v = Polynomial::new( vec![ 1.0 ] );
    let result = u.polydiv( &v );
    assert!( result.is_ok() );
    let ( q, r ) = result.unwrap();
    assert_eq!( q.degree().unwrap(), 2 );
    assert_eq!( q[0], 1.0 );
    assert_eq!( q[1], 2.0 );
    assert_eq!( q[2], 3.0 );
    assert_eq!( r.degree().unwrap(), 0 );
    assert_eq!( r[0], 0.0 );
}

#[test]
fn polydiv() {
    let u = Polynomial::new( vec![ 1.0, 1.0, 1.0 ] ); // 1 + x + x^2
    let v = Polynomial::new( vec![ 1.0, 1.0 ] ); // 1 + x
    let result = u.polydiv( &v );
    assert!( result.is_ok() );
    let ( q, r ) = result.unwrap();
    assert_eq!( q.degree().unwrap(), 1 );
    assert_eq!( q[0], 0.0 ); // q = x
    assert_eq!( q[1], 1.0 );
    assert_eq!( r.degree().unwrap(), 0 );
    assert_eq!( r[0], 1.0 ); // r = 1
}

#[test]
fn polydiv_cubic() {
    let u = Polynomial::cubic( 1.0, 0.0, 0.0, -1.0 ); // x^3 - 1
    let v = Polynomial::new( vec![ -1.0, 1.0 ] ); // -1 + x
    let result = u.polydiv( &v );
    assert!( result.is_ok() );
    let ( q, r ) = result.unwrap();
    assert_eq!( q.degree().unwrap(), 2 );
    assert_eq!( q[0], 1.0 ); // q = x^2 + x + 1
    assert_eq!( q[1], 1.0 );
    assert_eq!( q[2], 1.0 );
    assert_eq!( r.degree().unwrap(), 0 );
    assert_eq!( r[0], 0.0 ); // r = 0
}

#[test]
fn polynomial_addition_zero() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let zero = Polynomial::empty();
    let s = zero.clone() + p.clone();
    assert_eq!( s.degree().unwrap(), 2 );
    assert_eq!( s[0], 1.0 );
    assert_eq!( s[1], 2.0 );
    assert_eq!( s[2], 3.0 );
    let t = p + zero;
    assert_eq!( t.degree().unwrap(), 2 );
    assert_eq!( t[0], 1.0 );
    assert_eq!( t[1], 2.0 );
    assert_eq!( t[2], 3.0 );
}

#[test]
fn polynomial_addition() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let q = Polynomial::new( vec![ 4.0, 5.0, 6.0, 7.0 ] ); // 4 + 5x + 6x^2 + 7x^3
    let r = p.clone() + q;
    assert_eq!( r.degree().unwrap(), 3 );
    assert_eq!( r[0], 5.0 );
    assert_eq!( r[1], 7.0 );
    assert_eq!( r[2], 9.0 );
    assert_eq!( r[3], 7.0 );
    let zero = Polynomial::empty();
    let s = p + zero;
    assert_eq!( s.degree().unwrap(), 2 );
    assert_eq!( s[0], 1.0 );
    assert_eq!( s[1], 2.0 );
    assert_eq!( s[2], 3.0 );
}

#[test]
fn polynomial_negation() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let q = -p.clone();
    assert_eq!( q.degree().unwrap(), 2 );
    assert_eq!( q[0], -1.0 );
    assert_eq!( q[1], -2.0 );
    assert_eq!( q[2], -3.0 );
}

#[test]
fn polynomial_subtraction_zero() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let zero = Polynomial::empty();
    let s = zero.clone() - p.clone();
    assert_eq!( s.degree().unwrap(), 2 );
    assert_eq!( s[0], -1.0 );
    assert_eq!( s[1], -2.0 );
    assert_eq!( s[2], -3.0 );
    let t = p - zero;
    assert_eq!( t.degree().unwrap(), 2 );
    assert_eq!( t[0], 1.0 );
    assert_eq!( t[1], 2.0 );
    assert_eq!( t[2], 3.0 );
}

#[test]
fn polynomial_subtraction() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let q = Polynomial::new( vec![ 4.0, 5.0, 6.0, 7.0 ] ); // 4 + 5x + 6x^2 + 7x^3
    let r = p - q;
    assert_eq!( r.degree().unwrap(), 3 );
    assert_eq!( r[0], -3.0 );
    assert_eq!( r[1], -3.0 );
    assert_eq!( r[2], -3.0 );
    assert_eq!( r[3], -7.0 );
}

#[test]
fn polynomial_multiplication_zero() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let zero = Polynomial::empty();
    let s = zero.clone() * p.clone();
    match s.degree() {
        Ok( _ ) => assert!( false ),
        Err( _ ) => assert!( true )  // s has degree zero 
    }
    let t = p * zero;
    match t.degree() {
        Ok( _ ) => assert!( false ),
        Err( _ ) => assert!( true )  // t has degree zero 
    }
}

#[test]
fn polynomial_multiplication_one() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let one = Polynomial::new( vec![ 1.0 ] ); // 1
    let s = one.clone() * p.clone();
    assert_eq!( s.degree().unwrap(), 2 );
    assert_eq!( s[0], 1.0 );
    assert_eq!( s[1], 2.0 );
    assert_eq!( s[2], 3.0 );
    let t = p * one;
    assert_eq!( t.degree().unwrap(), 2 );
    assert_eq!( t[0], 1.0 );
    assert_eq!( t[1], 2.0 );
    assert_eq!( t[2], 3.0 );
}

#[test]
fn polynomial_multiplication() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2
    let q = Polynomial::new( vec![ 4.0, 5.0, 6.0, 7.0 ] ); // 4 + 5x + 6x^2 + 7x^3
    let r = p * q;
    assert_eq!( r.degree().unwrap(), 5 );
    assert_eq!( r[0], 4.0 );
    assert_eq!( r[1], 13.0 );
    assert_eq!( r[2], 28.0 );
    assert_eq!( r[3], 34.0 );
    assert_eq!( r[4], 32.0 );
    assert_eq!( r[5], 21.0 );
    let a = Polynomial::new( vec![ 5.0, 0.0, 10.0, 6.0 ] ); // 5 + 10x^2 + 6x^3
    let b = Polynomial::new( vec![ 1.0, 2.0, 4.0 ] ); // 1 + 2x + 4x^2
    let c = a * b;  
    assert_eq!( c.degree().unwrap(), 5 );
    assert_eq!( c[0], 5.0 );
    assert_eq!( c[1], 10.0 );
    assert_eq!( c[2], 30.0 );
    assert_eq!( c[3], 26.0 );
    assert_eq!( c[4], 52.0 );
    assert_eq!( c[5], 24.0 );
}

#[test]
fn polynomial_scalar_multiplication() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2 
    let s = p * 2.0;
    assert_eq!( s.degree().unwrap(), 2 );
    assert_eq!( s[0], 2.0 );
    assert_eq!( s[1], 4.0 );
    assert_eq!( s[2], 6.0 );
}

#[test]
fn remove_leading_zeros() {
    let mut p = Polynomial::new( vec![ 1.0, 2.0, 3.0, 0.0, 0.0 ] ); // 1 + 2x + 3x^2
    assert_eq!( p.degree().unwrap(), 4 );
    p.trim();
    assert_eq!( p.degree().unwrap(), 2 );
    assert_eq!( p[0], 1.0 );
    assert_eq!( p[1], 2.0 );
    assert_eq!( p[2], 3.0 );
}

#[test]
fn derivative() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2 
    let q = p.derivative();
    assert_eq!( q.degree().unwrap(), 1 );
    assert_eq!( q[0], 2.0 );
    assert_eq!( q[1], 6.0 );
    let r = q.derivative();
    assert_eq!( r.degree().unwrap(), 0 );
    assert_eq!( r[0], 6.0 );
}

#[test]
fn derivative_complex() { // 1 + (1+i)x + x^2 
    let p = Polynomial::new( vec![ Cmplx::one(), Cmplx::new( 1.0, 1.0 ),  Cmplx::one() ] ); 
    let q = p.derivative();
    assert_eq!( q.degree().unwrap(), 1 );
    assert_eq!( q[0], Cmplx::new( 1.0, 1.0 ) );
    assert_eq!( q[1], Cmplx::new( 2.0, 0.0 ) );
    let r = q.derivative();
    assert_eq!( r.degree().unwrap(), 0 );
    assert_eq!( r[0], Cmplx::new( 2.0, 0.0 ) );
}

#[test]
fn derivative_at() {
    let p = Polynomial::new( vec![ 1.0, 2.0, 3.0 ] ); // 1 + 2x + 3x^2 
    let q = p.derivative_at( 2.0, 1 );
    assert_eq!( q, 14.0 );
    let r = p.derivative_at( 2.0, 2 );
    assert_eq!( r, 6.0 );
}