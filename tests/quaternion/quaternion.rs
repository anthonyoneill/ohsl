use ohsl::quaternion::{Quaternion, Quat64};
use ohsl::traits::{Zero, One, Tiny};

#[test]
fn construction() {
    let q = Quaternion::<i32>::new( 1, 2, 3, 4 );
    assert_eq!( q.c[0], 1 );
    assert_eq!( q.c[1], 2 );
    assert_eq!( q.c[2], 3 );
    assert_eq!( q.c[3], 4 );
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    assert_eq!( q.c[0], 1.0 );
    assert_eq!( q.c[1], 2.0 );
    assert_eq!( q.c[2], 3.0 );
    assert_eq!( q.c[3], 4.0 );
}

#[test]
fn conjugate() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let qbar = q.conj();
    assert_eq!( qbar.c[0], 1.0 );
    assert_eq!( qbar.c[1], -2.0 );
    assert_eq!( qbar.c[2], -3.0 );
    assert_eq!( qbar.c[3], -4.0 );
}

#[test]
fn clone() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = q.clone();
    assert_eq!( r.c[0], 1.0 );
    assert_eq!( r.c[1], 2.0 );
    assert_eq!( r.c[2], 3.0 );
    assert_eq!( r.c[3], 4.0 );
}

#[test]
fn negation() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = -q;
    assert_eq!( r.c[0], -1.0 );
    assert_eq!( r.c[1], -2.0 );
    assert_eq!( r.c[2], -3.0 );
    assert_eq!( r.c[3], -4.0 );
}

#[test]
fn addition() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    let s = &q + &r; // Full non-consuming addition
    assert_eq!( s.c[0], 6.0 );
    assert_eq!( s.c[1], 8.0 );
    assert_eq!( s.c[2], 10.0 );
    assert_eq!( s.c[3], 12.0 );
    let t = q + &r; // Partial non-consuming addition
    assert_eq!( t.c[0], 6.0 );
    assert_eq!( t.c[1], 8.0 );
    assert_eq!( t.c[2], 10.0 );
    assert_eq!( t.c[3], 12.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let u = p + r;  // Consuming addition
    assert_eq!( u.c[0], 9.0 );
    assert_eq!( u.c[1], 9.0 );
    assert_eq!( u.c[2], 9.0 );
    assert_eq!( u.c[3], 9.0 );
}

#[test]
fn subtraction() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    let s = &q - &r; // Full non-consuming subtraction
    assert_eq!( s.c[0], -4.0 );
    assert_eq!( s.c[1], -4.0 );
    assert_eq!( s.c[2], -4.0 );
    assert_eq!( s.c[3], -4.0 );
    let t = q - &r; // Partial non-consuming subtraction
    assert_eq!( t.c[0], -4.0 );
    assert_eq!( t.c[1], -4.0 ); 
    assert_eq!( t.c[2], -4.0 );
    assert_eq!( t.c[3], -4.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let u = p - r;  // Consuming subtraction
    assert_eq!( u.c[0], -1.0 );
    assert_eq!( u.c[1], -3.0 );
    assert_eq!( u.c[2], -5.0 );
    assert_eq!( u.c[3], -7.0 );
}

#[test]
fn multiplication() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    let s = &q * &r; // Full non-consuming multiplication
    assert_eq!( s.c[0], -60.0 );
    assert_eq!( s.c[1], 12.0 );
    assert_eq!( s.c[2], 30.0 );
    assert_eq!( s.c[3], 24.0 );
    let w = &r * &q; // Full non-consuming multiplication (non-commutative)
    assert_eq!( w.c[0], -60.0 );
    assert_eq!( w.c[1], 20.0 );
    assert_eq!( w.c[2], 14.0 );
    assert_eq!( w.c[3], 32.0 );
    let t = q * &r; // Partial non-consuming multiplication
    assert_eq!( t.c[0], -60.0 );
    assert_eq!( t.c[1], 12.0 );
    assert_eq!( t.c[2], 30.0 );
    assert_eq!( t.c[3], 24.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let u = p * r;  // Consuming multiplication
    assert_eq!( u.c[0], -20.0 );
    assert_eq!( u.c[1], 48.0 );
    assert_eq!( u.c[2], 20.0 );
    assert_eq!( u.c[3], 46.0 );
}

#[test]
fn inverse() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let qinv = q.inv();
    assert_eq!( qinv.c[0], 0.03333333333333333 );
    assert_eq!( qinv.c[1], -0.06666666666666667 );
    assert_eq!( qinv.c[2], -0.10000000000000000 );
    assert_eq!( qinv.c[3], -0.13333333333333333 );
}

#[test]
fn division() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    let s = &q / &r; // Full non-consuming division
    let rinv = r.inv();
    let t = &q * &rinv; // Verify that q / r = q * r^-1
    assert!( (s.c[0] - t.c[0]).abs() < std::f64::EPSILON );
    assert!( (s.c[1] - t.c[1]).abs() < std::f64::EPSILON );
    assert!( (s.c[2] - t.c[2]).abs() < std::f64::EPSILON );
    assert!( (s.c[3] - t.c[3]).abs() < std::f64::EPSILON );
    let u = q / &r; // Partial non-consuming division
    assert!( (u.c[0] - t.c[0]).abs() < std::f64::EPSILON );
    assert!( (u.c[1] - t.c[1]).abs() < std::f64::EPSILON );
    assert!( (u.c[2] - t.c[2]).abs() < std::f64::EPSILON );
    assert!( (u.c[3] - t.c[3]).abs() < std::f64::EPSILON );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let w = p.clone() * rinv; // Verify that p / r = p * r^-1
    let v = p / r;  // Consuming division
    assert!( (v.c[0] - w.c[0]).abs() < std::f64::EPSILON );
    assert!( (v.c[1] - w.c[1]).abs() < std::f64::EPSILON );
    assert!( (v.c[2] - w.c[2]).abs() < std::f64::EPSILON );
    assert!( (v.c[3] - w.c[3]).abs() < std::f64::EPSILON );
}

#[test]
fn add_real() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = &q + &5.0; // Full non-consuming addition with a real number
    assert_eq!( r.c[0], 6.0 );
    assert_eq!( r.c[1], 2.0 );
    assert_eq!( r.c[2], 3.0 );
    assert_eq!( r.c[3], 4.0 );
    let s = q + &7.0; // Partial non-consuming addition with a real number
    assert_eq!( s.c[0], 8.0 );
    assert_eq!( s.c[1], 2.0 );
    assert_eq!( s.c[2], 3.0 );
    assert_eq!( s.c[3], 4.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let t = p + 11.0; // Consuming addition with a real number
    assert_eq!( t.c[0], 15.0 );
    assert_eq!( t.c[1], 3.0 );
    assert_eq!( t.c[2], 2.0 );
    assert_eq!( t.c[3], 1.0 );
}

#[test]
fn subtract_real() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = &q - &5.0; // Full non-consuming subtraction with a real number
    assert_eq!( r.c[0], -4.0 );
    assert_eq!( r.c[1], 2.0 );
    assert_eq!( r.c[2], 3.0 );
    assert_eq!( r.c[3], 4.0 );
    let s = q - &7.0; // Partial non-consuming subtraction with a real number
    assert_eq!( s.c[0], -6.0 );
    assert_eq!( s.c[1], 2.0 );
    assert_eq!( s.c[2], 3.0 );
    assert_eq!( s.c[3], 4.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let t = p - 11.0; // Consuming subtraction with a real number
    assert_eq!( t.c[0], -7.0 );
    assert_eq!( t.c[1], 3.0 );
    assert_eq!( t.c[2], 2.0 );
    assert_eq!( t.c[3], 1.0 );
}

#[test]
fn multiply_real() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = &q * &5.0; // Full non-consuming multiplication with a real number
    assert_eq!( r.c[0], 5.0 );
    assert_eq!( r.c[1], 10.0 );
    assert_eq!( r.c[2], 15.0 );
    assert_eq!( r.c[3], 20.0 );
    let s = q * &7.0; // Partial non-consuming multiplication with a real number
    assert_eq!( s.c[0], 7.0 );
    assert_eq!( s.c[1], 14.0 );
    assert_eq!( s.c[2], 21.0 );
    assert_eq!( s.c[3], 28.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let t = p * 11.0; // Consuming multiplication with a real number
    assert_eq!( t.c[0], 44.0 );
    assert_eq!( t.c[1], 33.0 );
    assert_eq!( t.c[2], 22.0 );
    assert_eq!( t.c[3], 11.0 );
    let u = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let v = 11.0 * u; // Multiplication on the left by a real number
    assert_eq!( v.c[0], 11.0 );
    assert_eq!( v.c[1], 22.0 );
    assert_eq!( v.c[2], 33.0 );
    assert_eq!( v.c[3], 44.0 );
}

#[test]
fn divide_real() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = &q / &5.0; // Full non-consuming division by a real number
    assert_eq!( r.c[0], 0.2 );
    assert_eq!( r.c[1], 0.4 );
    assert_eq!( r.c[2], 0.6 );
    assert_eq!( r.c[3], 0.8 );
    let s = q / &7.0; // Partial non-consuming division by a real number
    assert_eq!( s.c[0], 1.0 / 7.0 );
    assert_eq!( s.c[1], 2.0 / 7.0 );
    assert_eq!( s.c[2], 3.0 / 7.0 );
    assert_eq!( s.c[3], 4.0 / 7.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    let t = p / 11.0; // Consuming division by a real number
    assert_eq!( t.c[0], 4.0 / 11.0 );
    assert_eq!( t.c[1], 3.0 / 11.0 );
    assert_eq!( t.c[2], 2.0 / 11.0 );
    assert_eq!( t.c[3], 1.0 / 11.0 );
}

#[test]
fn add_assign() {
    let mut q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    q += &r; // Non-consuming addition assignment
    assert_eq!( q.c[0], 6.0 );
    assert_eq!( q.c[1], 8.0 );
    assert_eq!( q.c[2], 10.0 );
    assert_eq!( q.c[3], 12.0 );
    q += r; // Consuming addition assignment
    assert_eq!( q.c[0], 11.0 );
    assert_eq!( q.c[1], 14.0 );
    assert_eq!( q.c[2], 17.0 );
    assert_eq!( q.c[3], 20.0 );
    q += &5.0; // Non-consuming addition assignment with a real number
    assert_eq!( q.c[0], 16.0 );
    assert_eq!( q.c[1], 14.0 );
    assert_eq!( q.c[2], 17.0 );
    assert_eq!( q.c[3], 20.0 );
    q += 7.0; // Consuming addition assignment with a real number
    assert_eq!( q.c[0], 23.0 );
    assert_eq!( q.c[1], 14.0 );
    assert_eq!( q.c[2], 17.0 );
    assert_eq!( q.c[3], 20.0 );
}

#[test]
fn subtract_assign() {
    let mut q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    q -= &r; // Non-consuming subtraction assignment
    assert_eq!( q.c[0], -4.0 );
    assert_eq!( q.c[1], -4.0 );
    assert_eq!( q.c[2], -4.0 );
    assert_eq!( q.c[3], -4.0 );
    q -= r; // Consuming subtraction assignment
    assert_eq!( q.c[0], -9.0 );
    assert_eq!( q.c[1], -10.0 );
    assert_eq!( q.c[2], -11.0 );
    assert_eq!( q.c[3], -12.0 );
    q -= &5.0; // Non-consuming subtraction assignment with a real number
    assert_eq!( q.c[0], -14.0 );
    assert_eq!( q.c[1], -10.0 );
    assert_eq!( q.c[2], -11.0 );
    assert_eq!( q.c[3], -12.0 );
    q -= 7.0; // Consuming subtraction assignment with a real number
    assert_eq!( q.c[0], -21.0 );
    assert_eq!( q.c[1], -10.0 );
    assert_eq!( q.c[2], -11.0 );
    assert_eq!( q.c[3], -12.0 );
}

#[test]
fn multiply_assign() {
    let mut q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    q *= &r; // Non-consuming multiplication assignment
    assert_eq!( q.c[0], -60.0 );
    assert_eq!( q.c[1], 12.0 );
    assert_eq!( q.c[2], 30.0 );
    assert_eq!( q.c[3], 24.0 );
    let p = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    q *= p; // Consuming multiplication assignment
    assert_eq!( q.c[0], -360.0 );
    assert_eq!( q.c[1], -150.0 );
    assert_eq!( q.c[2], 60.0 );
    assert_eq!( q.c[3], -30.0 );
    q *= &11.0; // Non-consuming multiplication assignment with a real number
    assert_eq!( q.c[0], -3960.0 );
    assert_eq!( q.c[1], -1650.0 );
    assert_eq!( q.c[2], 660.0 );
    assert_eq!( q.c[3], -330.0 );
    q *= 7.0; // Consuming multiplication assignment with a real number
    assert_eq!( q.c[0], -27720.0 );
    assert_eq!( q.c[1], -11550.0 );
    assert_eq!( q.c[2], 4620.0 );
    assert_eq!( q.c[3], -2310.0 );
}

#[test]
fn divide_assign() {
    let mut q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = q.clone();
    let r = Quat64::new( 5.0, 6.0, 7.0, 8.0 );
    q /= &r; // Non-consuming division assignment
    let rinv = r.inv();
    let p = s * rinv; // Verify that q / r = q * r^-1
    assert!( (q.c[0] - p.c[0]).abs() < std::f64::EPSILON );
    assert!( (q.c[1] - p.c[1]).abs() < std::f64::EPSILON );
    assert!( (q.c[2] - p.c[2]).abs() < std::f64::EPSILON );
    assert!( (q.c[3] - p.c[3]).abs() < std::f64::EPSILON );
    let t = Quat64::new( 4.0, 3.0, 2.0, 1.0 );
    q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = q.clone();
    q /= t.clone(); // Consuming division assignment
    let u = s * t.inv(); // Verify that q / t = q * t^-1
    assert!( (q.c[0] - u.c[0]).abs() < std::f64::EPSILON );
    assert!( (q.c[1] - u.c[1]).abs() < std::f64::EPSILON );
    assert!( (q.c[2] - u.c[2]).abs() < std::f64::EPSILON );
    assert!( (q.c[3] - u.c[3]).abs() < std::f64::EPSILON );
    q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    q /= &11.0; // Non-consuming division assignment with a real number
    assert_eq!( q.c[0], 1.0 / 11.0 );
    assert_eq!( q.c[1], 2.0 / 11.0 );
    assert_eq!( q.c[2], 3.0 / 11.0 );
    assert_eq!( q.c[3], 4.0 / 11.0 );
    q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    q /= 7.0; // Consuming division assignment with a real number
    assert_eq!( q.c[0], 1.0 / 7.0 );
    assert_eq!( q.c[1], 2.0 / 7.0 );
    assert_eq!( q.c[2], 3.0 / 7.0 );
    assert_eq!( q.c[3], 4.0 / 7.0 );
}

#[test]
fn zero() {
    let q = Quat64::zero();
    assert_eq!( q.c[0], 0.0 );
    assert_eq!( q.c[1], 0.0 );
    assert_eq!( q.c[2], 0.0 );
    assert_eq!( q.c[3], 0.0 );
    let p = Quaternion::<i32>::zero();
    assert_eq!( p.c[0], 0 );
    assert_eq!( p.c[1], 0 );
    assert_eq!( p.c[2], 0 );
    assert_eq!( p.c[3], 0 );
    let r = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = &r + &q;
    assert_eq!( s.c[0], 1.0 );
    assert_eq!( s.c[1], 2.0 );
    assert_eq!( s.c[2], 3.0 );
    assert_eq!( s.c[3], 4.0 );
    let t = &r * &q;
    assert_eq!( t.c[0], 0.0 );
    assert_eq!( t.c[1], 0.0 );
    assert_eq!( t.c[2], 0.0 );
    assert_eq!( t.c[3], 0.0 );
}

#[test]
fn one() {
    let q = Quat64::one();
    assert_eq!( q.c[0], 1.0 );
    assert_eq!( q.c[1], 0.0 );
    assert_eq!( q.c[2], 0.0 );
    assert_eq!( q.c[3], 0.0 );
    let p = Quaternion::<i32>::one();
    assert_eq!( p.c[0], 1 );
    assert_eq!( p.c[1], 0 );
    assert_eq!( p.c[2], 0 );        
    assert_eq!( p.c[3], 0 );
    let r = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = &r + &q;
    assert_eq!( s.c[0], 2.0 );
    assert_eq!( s.c[1], 2.0 );
    assert_eq!( s.c[2], 3.0 );
    assert_eq!( s.c[3], 4.0 );
    let t = &q * &r;
    assert_eq!( t.c[0], 1.0 );
    assert_eq!( t.c[1], 2.0 );
    assert_eq!( t.c[2], 3.0 );
    assert_eq!( t.c[3], 4.0 );
    let u = &r * &q;
    assert_eq!( u.c[0], 1.0 );
    assert_eq!( u.c[1], 2.0 );
    assert_eq!( u.c[2], 3.0 );
    assert_eq!( u.c[3], 4.0 );
}

#[test]
fn tiny() {
    let q = Quat64::tiny();
    assert_eq!( q.c[0], 1.0e-40 );
    assert_eq!( q.c[1], 1.0e-40 );
    assert_eq!( q.c[2], 1.0e-40 );
    assert_eq!( q.c[3], 1.0e-40 );
}

#[test]
fn equality() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    assert!( q == r );
    let s = Quat64::new( 1.0, 2.0, 3.0, 4.0000000001 );
    assert!( q != s );
    let u = Quaternion::<i32>::new( 1, 2, 3, 4 );
    let v = Quaternion::<i32>::new( 1, 2, 3, 4 );
    assert!( u == v );
    let w = Quaternion::<i32>::new( 1, 2, 7, 4 );
    assert!( u != w );  
} 

#[test]
fn ordering() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    assert!( q == r );
    let s = Quat64::new( 1.0, 2.0, 3.0, 5.0 );
    assert!( s > q );
    let t = Quat64::new( 1.0, 2.0, 3.0, 3.0 );
    assert!( t < q );
    let u = Quat64::new( 1.0, 2.0, 4.0, 0.0 );
    assert!( u > q );
    let v = Quat64::new( 1.0, 2.0, 2.0, 5.0 );
    assert!( v < q );
    let w = Quat64::new( 1.0, 3.0, 0.0, 0.0 );
    assert!( w > q );
    let x = Quat64::new( 0.0, 3.0, 0.0, 0.0 );
    assert!( x < q );
    let y = Quat64::new( 2.0, 0.0, 0.0, 0.0 );
    assert!( y > q );
    let z = Quat64::new( 0.0, 0.0, 0.0, 0.0 );
    assert!( z < q );
    let a = Quat64::new( -1.0, 3.0, 4.0, 5.0 );
    assert!( a < q );
    let b = Quat64::new( 1.0, 2.0, -1.0, 4.0 );
    assert!( b < q );
    let c = Quat64::new( 1.0, 2.0, 3.0, -1.0 );
    assert!( c < q );
}

#[test]
fn abs_sqr() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = q.abs_sqr();
    assert_eq!( r, 30.0 );
    let p = Quaternion::<i32>::new( 5, 7, 9, 11 );
    let s = p.abs_sqr();
    assert_eq!( s, 276 );
}

#[test]
fn norm() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let r = q.norm();
    assert_eq!( r, 30.0_f64.sqrt() );
    let p = Quat64::new( 5.0, 7.0, 9.0, 11.0 );
    let s = p.norm();
    assert_eq!( s, 276.0_f64.sqrt() );
}

#[test]
fn display() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = format!("{}", q);
    assert_eq!( s, "1 + 2i + 3j + 4k" );
    let p = Quaternion::<i32>::new( 5, 7, 9, 11 );
    let t = format!("{}", p);
    assert_eq!( t, "5 + 7i + 9j + 11k" );
}

#[test]
fn debug() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = format!("{:?}", q);
    assert_eq!( s, "Quaternion { c: [1.0, 2.0, 3.0, 4.0] }" );
    let p = Quaternion::<i32>::new( 5, 7, 9, 11 );
    let t = format!("{:?}", p);
    assert_eq!( t, "Quaternion { c: [5, 7, 9, 11] }" );
}  

#[test]
fn scalar_and_vector_parts() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let s = q.scalar();
    assert_eq!( s, 1.0 );
    let v = q.vector();
    assert_eq!( v[0], 2.0 );
    assert_eq!( v[1], 3.0 );
    assert_eq!( v[2], 4.0 );
    assert_eq!( v.len(), 3 );
    let p = Quaternion::<i32>::new( 5, 7, 9, 11 );
    let t = p.scalar();
    assert_eq!( t, 5 );
    let w = p.vector();
    assert_eq!( w, [7, 9, 11] );
}

#[test]
fn versor() {
    let q = Quat64::new( 1.0, 2.0, 3.0, 4.0 );
    let v = q.versor();
    assert_eq!( v.c[0], 1.0 / 30.0_f64.sqrt() );
    assert_eq!( v.c[1], 2.0 / 30.0_f64.sqrt() );
    assert_eq!( v.c[2], 3.0 / 30.0_f64.sqrt() );
    assert_eq!( v.c[3], 4.0 / 30.0_f64.sqrt() );
    let norm = v.norm();
    assert!((norm - 1.0).abs() < std::f64::EPSILON );

}