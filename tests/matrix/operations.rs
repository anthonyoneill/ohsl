use ohsl::matrix::Matrix;
use ohsl::vector::Vector;

fn m_matrix() -> Matrix::<i32> {
    let mut m = Matrix::<i32>::new( 2, 3, 1 );
    m[(0,1)] = 2;
    m[(0,2)] = 3; 
    m[(1,0)] = 4;
    m[(1,1)] = 5;
    m[(1,2)] = 6;
    m
}

#[test]
fn get_row() {
    let m = m_matrix();
    let v: Vector<i32> = m.get_row( 0 );
    assert_eq!( v[0], 1 );
    assert_eq!( v[1], 2 );
    assert_eq!( v[2], 3 );
}

#[test]
fn get_col() {
    let m = m_matrix();
    let u: Vector<i32> = m.get_col( 0 );
    assert_eq!( u[0], 1 );
    assert_eq!( u[1], 4 );
}

#[test]
fn set_row() {
    let mut m = m_matrix();
    let new_row = Vector::<i32>::new( 3, 7 );
    m.set_row( 0, new_row );
    assert_eq!( m[(0,0)], 7 );
    assert_eq!( m[(0,1)], 7 );
    assert_eq!( m[(0,2)], 7 );
}

#[test]
fn set_col() {
    let mut m = m_matrix();
    let new_col = Vector::<i32>::new( 2, 8 );
    m.set_col( 1, new_col );
    assert_eq!( m[(0,1)], 8 );
    assert_eq!( m[(1,1)], 8 );
}

#[test]
fn delete_row() {
    let mut m = m_matrix();
    m.delete_row( 0 );
    assert_eq!( m[(0,0)], 4 );
    assert_eq!( m[(0,1)], 5 );
    assert_eq!( m[(0,2)], 6 );
    assert_eq!( m.rows(), 1 );
}

#[test]
fn resize() {
    let mut m = Matrix::<i32>::new( 2, 3, 1 );
    m[(0,1)] = 2;
    m[(0,2)] = 3; 
    m[(1,0)] = 4;
    m[(1,1)] = 5;
    m[(1,2)] = 6;
    m.resize( 2, 2 );
    assert_eq!( m.rows(), 2 );
    assert_eq!( m.cols(), 2 );
    assert_eq!( m[(0,0)], 1 );
    m.resize( 1, 4 );
    assert_eq!( m.rows(), 1 );
    assert_eq!( m.cols(), 4 );
    assert_eq!( m[(0,1)], 2 );
    assert_eq!( m[(0,3)], 0 );
}

#[test]
fn transpose_in_place() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[(0,1)] = 2;
    m[(1,0)] = 3; 
    m[(1,1)] = 4;
    m.transpose_in_place();
    assert_eq!( m[(0,1)], 3 );
    assert_eq!( m[(1,0)], 2 );
    let mut n = Matrix::<i32>::new( 3, 2, 1 );
    n[(0,1)] = 2;
    n[(1,0)] = 3;
    n[(1,1)] = 4;
    n[(2,0)] = 5;
    n[(2,1)] = 6;
    n.transpose_in_place();
    assert_eq!( n.rows(), 2 );
    assert_eq!( n.cols(), 3 );
    assert_eq!( n[(0,1)], 3 );
    assert_eq!( n[(0,2)], 5 );
    assert_eq!( n[(1,2)], 6 );
}

#[test]
fn test_matrix_transpose() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[(0,1)] = 2;
    m[(1,0)] = 3; 
    m[(1,1)] = 4;
    let p = m.transpose();
    assert_eq!( p[(0,1)], 3 );
    assert_eq!( p[(1,0)], 2 );
}

#[test]
fn fill() {
    let mut m = Matrix::<i32>::new( 3, 3, 1 );
    m.fill( 0 );
    assert_eq!( m[(0,0)], 0 );
    let mut n = Matrix::<i32>::new( 4, 3, 0 );
    n.fill_diag( 3 );
    assert_eq!( n[(2,2)], 3 );
    assert_eq!( n[(3,2)], 0 );
    n.fill_band( -1, 7 );
    assert_eq!( n[(2,1)], 7 );
    assert_eq!( n[(3,2)], 7 );
    n.fill_tridiag( 1, 2, 3 );
    assert_eq!( n[(0,0)], 2 );
    assert_eq!( n[(1,0)], 1 );
    assert_eq!( n[(0,1)], 3 );
}

#[test]
fn swap() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    m[(0,1)] = 2;
    m[(1,0)] = 3; 
    m[(1,1)] = 4;
    m.swap_rows( 0, 1 );
    assert_eq!( m[(0,0)], 3 );
    assert_eq!( m[(1,1)], 2 );
    m.swap_elem( 0, 0, 1, 1 );
    assert_eq!( m[(0,0)], 2 );
    assert_eq!( m[(1,1)], 3 );
}

