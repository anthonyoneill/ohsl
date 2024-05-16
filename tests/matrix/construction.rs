use ohsl::matrix::Matrix;

#[test]
fn unspecified_size() {
    let m = Matrix::<i32>::empty();
    assert_eq!( m.rows(), 0 );
    assert_eq!( m.cols(), 0 );
}

#[test]
fn specified_size() {
    let mut m = Matrix::<i32>::new( 2, 2, 1 );
    assert_eq!( m.rows(), 2 );
    assert_eq!( m.cols(), 2 );
    assert_eq!( m[(0,0)], 1 );
    assert_eq!( m[(0,1)], 1 );
    assert_eq!( m[(1,0)], 1 );
    assert_eq!( m[(1,1)], 1 );
    m[(0,0)] = 7;
    assert_eq!( m[(0,0)], 7 );
}

#[test]
fn clone() {
    let m = Matrix::<i32>::new( 3, 3, 7 );
    let n = m.clone();
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!( m[(i,j)], n[(i,j)] );
        }
    };
}

#[test]
fn test_matrix_eye() {
    let eye = Matrix::<i32>::eye( 3 );
    assert_eq!( eye.rows(), 3 );
    assert_eq!( eye.cols(), 3 );
    assert_eq!( eye[(0,1)], 0 );
    assert_eq!( eye[(2,2)], 1 );
}