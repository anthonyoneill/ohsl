use ohsl::matrix::Matrix;
use ohsl::vector::Vector;

#[test]
fn unary_minus() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    assert_eq!( m[1][1], 1 );
    let w = -m; 
    assert_eq!( w[1][1], -1 );
}

#[test]
fn binary_plus() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    let n = Matrix::<i32>::new( 2, 2, 4 );
    let p = m + n; 
    assert_eq!( p[1][1], 5 );
}

#[test]
fn binary_minus() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    let n = Matrix::<i32>::new( 2, 2, 4 );
    let q = m - n;
    assert_eq!( q[1][1], -3 );
}

#[test]
fn scalar_multiply() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    let n = m * 6;
    assert_eq!( n[0][1], 6 );
    let m = Matrix::<f64>::new( 2, 2, 1.4 );
    let n = 2.0 * m;
    assert_eq!( n[0][0], 2.8 );
}

#[test]
fn scalar_divide() {
    let m = Matrix::<f64>::new( 2, 2, 1.4 );
    let p = m / 2.0;
    assert_eq!( p[1][1], 0.7 );
}

#[test]
fn add_assign() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    let v = Matrix::<i32>::new( 2, 2, 3 );
    m += v;
    assert_eq!( m[0][0], 4 );
}

#[test]
fn subtract_assign() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    let u = Matrix::<i32>::new( 2, 2, 5 );
    m -= u;
    assert_eq!( m[1][1], -4 );
}

#[test]
fn scalar_assign() {
    let mut m = Matrix::<f64>::new( 2, 2, 1.0 );
    m *= 4.0;
    assert_eq!( m[0][0], 4.0 );
    m /= 8.0;
    assert_eq!( m[1][1], 0.5 );
    m += 1.0;
    assert_eq!( m[1][0], 1.5 );
    m -= 0.5;
    assert_eq!( m[0][1], 1.0 );
}

#[test]
fn matrix_multiplication() {
    let mut m = Matrix::<i32>::new( 2, 3, 1 );
    m[0][1] = 2;
    m[0][2] = 3; 
    m[1][0] = 4;
    m[1][1] = 5;
    m[1][2] = 6;
    let mut n = Matrix::<i32>::new( 3, 2, 1 );
    n[0][0] = 5;
    n[0][1] = 6; 
    n[1][0] = 4;
    n[1][1] = 3;
    n[2][0] = 2;
    n[2][1] = 1;
    let v = Vector::<i32>::new( 3, 1 );
    let u: Vector<i32> = &m * &v; // Non-consuming matrix-vector multiplication
    assert_eq!( u[0], 6 );
    assert_eq!( u[1], 15 );
    let p: Matrix<i32> = &m * &n; // Non-consuming matrix-matrix multiplication
    assert_eq!( p[0][0], 19 );
    assert_eq!( p[0][1], 15 );
    assert_eq!( p[1][0], 52 );
    assert_eq!( p[1][1], 45 );
    m.clear();
    assert_eq!( m.rows(), 0 );
    assert_eq!( m.cols(), 0 );
}