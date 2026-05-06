use ohsl::triad::Triad;
use ohsl::vector::Vec64;

#[test]
fn unspecified_size() {
    let t = Triad::<i32>::empty();
    assert_eq!( t.panels(), 0 );
    assert_eq!( t.rows(), 0 );
    assert_eq!( t.cols(), 0 );
}

#[test]
fn specified_size() {
    let mut t = Triad::<i32>::new( 2, 2, 2, 1 );
    assert_eq!( t.panels(), 2 );
    assert_eq!( t.rows(), 2 );
    assert_eq!( t.cols(), 2 );
    assert_eq!( t[(0,0,0)], 1 );
    assert_eq!( t[(0,0,1)], 1 );
    assert_eq!( t[(0,1,0)], 1 );
    assert_eq!( t[(0,1,1)], 1 );
    assert_eq!( t[(1,0,0)], 1 );
    assert_eq!( t[(1,0,1)], 1 );
    assert_eq!( t[(1,1,0)], 1 );
    assert_eq!( t[(1,1,1)], 1 );
    t[(0,0,0)] = 7;
    assert_eq!( t[(0,0,0)], 7 );
}

#[test]
fn from_vec() {
    let vec = vec![ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ];
    let t = Triad::<i32>::from_vec( 2, 3, 2, vec.clone() );
    assert_eq!( t.panels(), 2 );
    assert_eq!( t.rows(), 3 );
    assert_eq!( t.cols(), 2 );
    assert_eq!( t[(0,0,0)], 1 );
    assert_eq!( t[(0,0,1)], 2 );
    assert_eq!( t[(0,1,0)], 3 );
    assert_eq!( t[(0,1,1)], 4 );
    assert_eq!( t[(0,2,0)], 5 );
    assert_eq!( t[(0,2,1)], 6 );
    
    assert_eq!( t[(1,0,0)], 7 );
    assert_eq!( t[(1,0,1)], 8 );
    assert_eq!( t[(1,1,0)], 9 );
    assert_eq!( t[(1,1,1)], 10 );
    assert_eq!( t[(1,2,0)], 11 );
    assert_eq!( t[(1,2,1)], 12 );
    let u = Triad::<i32>::from_vec( 3, 2, 2, vec );
    assert_eq!( u.panels(), 3 );
    assert_eq!( u.rows(), 2 );
    assert_eq!( u.cols(), 2 );
    assert_eq!( u[(0,0,0)], 1 );
    assert_eq!( u[(0,0,1)], 2 );
    assert_eq!( u[(0,1,0)], 3 );
    assert_eq!( u[(0,1,1)], 4 );

    assert_eq!( u[(1,0,0)], 5 );
    assert_eq!( u[(1,0,1)], 6 );
    assert_eq!( u[(1,1,0)], 7 );
    assert_eq!( u[(1,1,1)], 8 );

    assert_eq!( u[(2,0,0)], 9 );
    assert_eq!( u[(2,0,1)], 10 );
    assert_eq!( u[(2,1,0)], 11 );
    assert_eq!( u[(2,1,1)], 12 );
    let vec = vec![ 
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 
    ];
    let v = Triad::<i32>::from_vec( 4, 3, 2, vec );
    assert_eq!( v.panels(), 4 );
    assert_eq!( v.rows(), 3 );
    assert_eq!( v.cols(), 2 );
    assert_eq!( v[(0,0,0)], 1 );
    assert_eq!( v[(0,0,1)], 2 );
    assert_eq!( v[(0,1,0)], 3 );
    assert_eq!( v[(0,1,1)], 4 );
    assert_eq!( v[(0,2,0)], 5 );
    assert_eq!( v[(0,2,1)], 6 );

    assert_eq!( v[(1,0,0)], 7 );
    assert_eq!( v[(1,0,1)], 8 );
    assert_eq!( v[(1,1,0)], 9 );
    assert_eq!( v[(1,1,1)], 10 );
    assert_eq!( v[(1,2,0)], 11 );
    assert_eq!( v[(1,2,1)], 12 );

    assert_eq!( v[(2,0,0)], 13 );
    assert_eq!( v[(2,0,1)], 14 );
    assert_eq!( v[(2,1,0)], 15 );
    assert_eq!( v[(2,1,1)], 16 );
    assert_eq!( v[(2,2,0)], 17 );
    assert_eq!( v[(2,2,1)], 18 );

    assert_eq!( v[(3,0,0)], 19 );
    assert_eq!( v[(3,0,1)], 20 );
    assert_eq!( v[(3,1,0)], 21 );
    assert_eq!( v[(3,1,1)], 22 );
    assert_eq!( v[(3,2,0)], 23 );
    assert_eq!( v[(3,2,1)], 24 );
}

#[test]
fn clone() {
    let t = Triad::<i32>::new( 2, 2, 2, 7 );
    let u = t.clone();
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                assert_eq!( t[(i,j,k)], u[(i,j,k)] );
            }
        }
    };
}

fn vector_function( x: &Vec64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = f64::powf( x[0], 3.0 ) + x[1] -1.0;
    f[1] = f64::powf( x[1], 3.0 ) - x[0] + 1.0;
    f
}

#[test]
fn hessian() {
    let point = Vec64::create( vec![ 1.0, 1.0 ] );
    let hes = Triad::<f64>::hessian( &point, &vector_function, 1e-6 );
    assert_eq!( hes.panels(), 2 );
    assert_eq!( hes.rows(), 2 );
    assert_eq!( hes.cols(), 2 );
    assert!( ( hes[(0,0,0)] - 6.0 ).abs() < 1e-4 );
    assert!( hes[(0,0,1)].abs() < 1e-4 );
    assert!( hes[(0,1,0)].abs() < 1e-4 );
    assert!( hes[(0,1,1)].abs() < 1e-4 );
    assert!( hes[(1,0,0)].abs() < 1e-4 );
    assert!( hes[(1,0,1)].abs() < 1e-4 );
    assert!( hes[(1,1,0)].abs() < 1e-4 );
    assert!( (hes[(1,1,1)] - 6.0 ).abs() < 1e-4 );
}