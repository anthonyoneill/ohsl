use ohsl::matrix::{Matrix, Mat64};
use ohsl::vector::{Vector, Vec64};
use ohsl::complex::Cmplx;

#[test]
fn test_matrix_unspecified_size() {
    let m = Matrix::<i32>::empty();
    assert_eq!( m.rows(), 0 );
    assert_eq!( m.cols(), 0 );
}

#[test]
fn test_matrix_specified_size() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    assert_eq!( m.rows(), 2 );
    assert_eq!( m.cols(), 2 );
    assert_eq!( m[0][0], 1 );
    assert_eq!( m[0][1], 1 );
    assert_eq!( m[1][0], 1 );
    assert_eq!( m[1][1], 1 );
    m[0][0] = 7;
    assert_eq!( m[0][0], 7 );
}

#[test]
fn test_create_clone() {
    let v1 = Vector::<i32>::new( 3, 3 );
    let v2 = Vector::<i32>::new( 3, 7 );
    let mut vec = Vec::new();
    vec.push( v1 );
    vec.push( v2 );
    let m = Matrix::<i32>::create( vec );
    assert_eq!( m[0][0], 3 );
    assert_eq!( m[1][2], 7 );
    assert_eq!( m.rows(), 2 );
    assert_eq!( m.cols(), 3 );
    assert_eq!( m.numel(), 6 );
    let n = m.clone();
    assert_eq!( n[1][2], 7 );
}

#[test]
fn test_matrix_unary_minus() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    assert_eq!( m[1][1], 1 );
    let w = -m; 
    assert_eq!( w[1][1], -1 );
}

#[test]
fn test_matrix_binary_plus_minus() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    let n = Matrix::<i32>::new( 2, 2, 4 );
    let p = m.clone() + n.clone(); 
    assert_eq!( p[1][1], 5 );
    let q = m - n;
    assert_eq!( q[1][1], -3 )
}

#[test]
fn test_matrix_scalar_mul_div() {
    let m = Matrix::<i32>::new( 2, 2, 1 );
    let n = m * 6;
    assert_eq!( n[0][1], 6 );
    let m = Matrix::<f64>::new( 2, 2, 1.4 );
    let n = 2.0 * m;
    assert_eq!( n[0][0], 2.8 );
    let p = n / 4.0;
    assert_eq!( p[1][1], 0.7 );
}

#[test]
fn test_matrix_add_sub_assign() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    let v = Matrix::<i32>::new( 2, 2, 3 );
    m += v;
    assert_eq!( m[0][0], 4 );
    let u = Matrix::<i32>::new( 2, 2, 5 );
    m -= u;
    assert_eq!( m[1][1], -1 );
}

#[test]
fn test_matrix_scalar_assign() {
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
fn test_matrix_get_set() {
    let mut m = Matrix::<i32>::new( 2, 3, 1 );
    m[0][1] = 2;
    m[0][2] = 3; 
    m[1][0] = 4;
    m[1][1] = 5;
    m[1][2] = 6;
    let v: Vector<i32> = m.get_row( 0 );
    assert_eq!( v[0], 1 );
    assert_eq!( v[1], 2 );
    assert_eq!( v[2], 3 );
    let u: Vector<i32> = m.get_col( 0 );
    assert_eq!( u[0], 1 );
    assert_eq!( u[1], 4 );
    let new_row = Vector::<i32>::new( 3, 7 );
    m.set_row( 0, new_row );
    assert_eq!( m[0][0], 7 );
    assert_eq!( m[0][1], 7 );
    assert_eq!( m[0][2], 7 );
    let new_col = Vector::<i32>::new( 2, 8 );
    m.set_col( 1, new_col );
    assert_eq!( m[0][1], 8 );
    assert_eq!( m[1][1], 8 );
    m.delete_row( 0 );
    assert_eq!( m[0][0], 4 );
    assert_eq!( m[0][1], 8 );
    assert_eq!( m[0][2], 6 );
    assert_eq!( m.rows(), 1 );
}

#[test]
fn test_matrix_multiplication() {
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
    let u: Vector<i32> = m.clone() * v;
    assert_eq!( u[0], 6 );
    assert_eq!( u[1], 15 );
    let p = m.clone() * n;
    assert_eq!( p[0][0], 19 );
    assert_eq!( p[0][1], 15 );
    assert_eq!( p[1][0], 52 );
    assert_eq!( p[1][1], 45 );
    m.clear();
    assert_eq!( m.rows(), 0 );
    assert_eq!( m.cols(), 0 );
}

#[test]
fn test_matrix_eye() {
    let eye = Matrix::<i32>::eye( 3 );
    assert_eq!( eye.rows(), 3 );
    assert_eq!( eye.cols(), 3 );
    assert_eq!( eye[0][1], 0 );
    assert_eq!( eye[2][2], 1 );
}

#[test]
fn test_matrix_resize() {
    let mut m = Matrix::<i32>::new( 2, 3, 1 );
    m[0][1] = 2;
    m[0][2] = 3; 
    m[1][0] = 4;
    m[1][1] = 5;
    m[1][2] = 6;
    m.resize( 2, 2 );
    assert_eq!( m.rows(), 2 );
    assert_eq!( m.cols(), 2 );
    assert_eq!( m[0][0], 1 );
    m.resize( 1, 4 );
    assert_eq!( m.rows(), 1 );
    assert_eq!( m.cols(), 4 );
    assert_eq!( m[0][1], 2 );
    assert_eq!( m[0][3], 0 );
}

#[test]
fn test_matrix_transpose() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[0][1] = 2;
    m[1][0] = 3; 
    m[1][1] = 4;
    m.transpose_in_place();
    assert_eq!( m[0][1], 3 );
    assert_eq!( m[1][0], 2 );
    let mut n = Matrix::<i32>::new( 3, 2, 1 );
    n[0][1] = 2;
    n[1][0] = 3;
    n[1][1] = 4;
    n[2][0] = 5;
    n[2][1] = 6;
    n.transpose_in_place();
    assert_eq!( n.rows(), 2 );
    assert_eq!( n.cols(), 3 );
    assert_eq!( n[0][1], 3 );
    assert_eq!( n[0][2], 5 );
    assert_eq!( n[1][2], 6 );
    let p = m.transpose();
    assert_eq!( p[0][1], 2 );
    assert_eq!( p[1][0], 3 );
}

#[test]
fn test_matrix_fill() {
    let mut m = Matrix::<i32>::new( 3, 3, 1 );
    m.fill( 0 );
    assert_eq!( m[0][0], 0 );
    let mut n = Matrix::<i32>::new( 4, 3, 0 );
    n.fill_diag( 3 );
    assert_eq!( n[2][2], 3 );
    assert_eq!( n[3][2], 0 );
    n.fill_band( -1, 7 );
    assert_eq!( n[2][1], 7 );
    assert_eq!( n[3][2], 7 );
    n.fill_tridiag( 1, 2, 3 );
    assert_eq!( n[0][0], 2 );
    assert_eq!( n[1][0], 1 );
    assert_eq!( n[0][1], 3 );
}

#[test]
fn test_matrix_swap() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[0][1] = 2;
    m[1][0] = 3; 
    m[1][1] = 4;
    m.swap_rows( 0, 1 );
    assert_eq!( m[0][0], 3 );
    assert_eq!( m[1][1], 2 );
    m.swap_elem( 0, 0, 1, 1 );
    assert_eq!( m[0][0], 2 );
    assert_eq!( m[1][1], 3 );
}

#[test]
fn test_matrix_norms() {
    let mut a = Mat64::new( 2, 2, 1.0 );
    a[0][1] = -7.0;
    a[1][0] = -2.0;
    a[1][1] = -3.0;
    assert_eq!( a.norm_1(), 10.0 );
    assert_eq!( a.norm_inf(), 8.0 );
    let rounded_to_3dp = (a.norm_p( 2.0 ) * 1000.0).round() / 1000.0;
    assert_eq!( rounded_to_3dp, 7.937 );
    let difference = a.norm_frob() - 7.937;
    assert!( difference.abs() < 1e-3 );
    assert_eq!( a.norm_max(), 7.0 );
}

/*#[test]
fn test_matrix_file_output() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[0][1] = 2;
    m[1][0] = 3; 
    m[1][1] = 4;
    m.output( "./output.txt" );
}*/

#[test]
fn test_matrix_solve_basic() {
    let mut a = Mat64::new( 2, 2, 1.0 );
    a[0][1] = 2.0;
    a[1][0] = 3.0; 
    a[1][1] = 4.0;
    let b = Vec64::create( vec![ 5.0, 11.0 ] ); 
    let x = a.solve_basic( b );
    assert_eq!( x[0], 1.0 );
    assert_eq!( x[1], 2.0 );
    let mut c = Matrix::<Cmplx>::new( 2, 2, Cmplx::new( 1.0, 1.0 ) );
    c[0][1] = Cmplx::new( -1.0, 0.0 );
    c[1][0] = Cmplx::new( 1.0, -1.0 );
    let mut d = Vector::<Cmplx>::new( 2, Cmplx::new( 0.0, 1.0 ) );
    d[1] = Cmplx::new( 1.0, 0.0 );
    let y = c.solve_basic( d );
    assert_eq!( y[0], Cmplx::new( 0.5, 0.5 ) );
    assert_eq!( y[1], Cmplx::new( 0.0, 0.0 ) );
}

#[test]
fn test_matrix_solve_lu() {
    let mut a = Mat64::new( 2, 2, 4.0 );
    a[0][1] = 3.0;
    a[1][0] = 6.0; 
    a[1][1] = 3.0;
    //let (pivots, P) = a.lu_decomp_in_place();
    let b = Vec64::create( vec![ 10.0, 12.0 ] );
    let x = a.solve_lu( b );
    assert_eq!( x[0], 1.0 );
    assert_eq!( x[1], 2.0 );
}

#[test]
fn test_matrix_determinant() {
    let mut a = Mat64::new( 2, 2, 4.0 );
    a[0][1] = 3.0;
    a[1][0] = 6.0; 
    a[1][1] = 3.0;
    let det = a.determinant();
    assert_eq!( det, -6.0 );
}

#[test]
fn test_matrix_inverse() {
    let mut a = Mat64::new( 2, 2, 4.0 );
    a[0][1] = 3.0;
    a[1][0] = 6.0; 
    a[1][1] = 3.0;
    let inv = a.inverse();
    assert_eq!( inv[0][0], -0.5 );
    assert_eq!( inv[0][1], 0.5 );
    assert_eq!( inv[1][0], 1.0 );
    assert_eq!( inv[1][1], -2.0/3.0 );
}