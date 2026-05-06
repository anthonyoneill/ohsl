use ohsl::triad::Triad;
use ohsl::vector::Vector;

#[test]
fn get_panel() {
    let vec = vec![ 
        1, 2, 
        3, 4, 
        
        5, 6, 
        7, 8, 
        
        9, 10, 
        11, 12 
    ];
    let u = Triad::<i32>::from_vec( 3, 2, 2, vec );
    assert_eq!( u.panels(), 3 );
    assert_eq!( u.rows(), 2 );
    assert_eq!( u.cols(), 2 );
    let panel0 = u.get_panel( 0 );
    assert_eq!( panel0.rows(), 2 );
    assert_eq!( panel0.cols(), 2 );
    assert_eq!( panel0[(0,0)], 1 );
    assert_eq!( panel0[(0,1)], 2 );
    assert_eq!( panel0[(1,0)], 3 );
    assert_eq!( panel0[(1,1)], 4 );
    let panel1 = u.get_panel( 1 );
    assert_eq!( panel1.rows(), 2 );
    assert_eq!( panel1.cols(), 2 );
    assert_eq!( panel1[(0,0)], 5 );
    assert_eq!( panel1[(0,1)], 6 );
    assert_eq!( panel1[(1,0)], 7 );
    assert_eq!( panel1[(1,1)], 8 );
    let panel2 = u.get_panel( 2 );
    assert_eq!( panel2.rows(), 2 );
    assert_eq!( panel2.cols(), 2 );
    assert_eq!( panel2[(0,0)], 9 );
    assert_eq!( panel2[(0,1)], 10 );
    assert_eq!( panel2[(1,0)], 11 );
    assert_eq!( panel2[(1,1)], 12 );
}

#[test]
fn multiply() {
    let x = Vector::<i32>::create( vec![ 1, 2, 3 ] );
    let b = Triad::<i32>::from_vec( 3, 2, 3, vec![
        1, 2, 3,
        4, 5, 6,

        7, 8, 9,
        10, 11, 12,

        13, 14, 15,
        16, 17, 18
    ] );
    let a = b.multiply( &x );
    assert_eq!( a.rows(), 3 );
    assert_eq!( a.cols(), 2 );
    assert_eq!( a[(0,0)], 14 );
    assert_eq!( a[(0,1)], 32 );
    assert_eq!( a[(1,0)], 50 );
    assert_eq!( a[(1,1)], 68 );
    assert_eq!( a[(2,0)], 86 );
    assert_eq!( a[(2,1)], 104 );

    let y = Vector::<i32>::create( vec![ 1, 2 ] );
    let c = Triad::<i32>::from_vec( 2, 3, 2, vec![
        1, 2,
        3, 4,
        5, 6,

        7, 8,
        9, 10,
        11, 12
    ] );
    let d = c.multiply( &y );
    assert_eq!( d.rows(), 2 );
    assert_eq!( d.cols(), 3 );
    assert_eq!( d[(0,0)], 5 );
    assert_eq!( d[(0,1)], 11 );
    assert_eq!( d[(0,2)], 17 );
    assert_eq!( d[(1,0)], 23 );
    assert_eq!( d[(1,1)], 29 );
    assert_eq!( d[(1,2)], 35 );
}

#[test]
fn multiply_twice() {
    let x = Vector::<i32>::create( vec![ 1, 2 ] );
    let b = Triad::<i32>::from_vec( 2, 2, 2, vec![
        1, 2,
        3, 4,

        5, 6,
        7, 8
    ] );
    let a = b.multiply( &x );
    assert_eq!( a.rows(), 2 );
    assert_eq!( a.cols(), 2 );
    assert_eq!( a[(0,0)], 5 );
    assert_eq!( a[(0,1)], 11 );
    assert_eq!( a[(1,0)], 17 );
    assert_eq!( a[(1,1)], 23 );
    let c = &a * &x;
    assert_eq!( c.size(), 2 );
    assert_eq!( c[0], 27 );
    assert_eq!( c[1], 63 );
    let d = b.multiply_twice( &x );
    assert_eq!( d.size(), 2 );
    assert_eq!( d[0], 27 );
    assert_eq!( d[1], 63 );

    let y = Vector::<i32>::create( vec![ 1, 2, 3 ] );
    let e = Triad::<i32>::from_vec( 3, 3, 3, vec![
        1, 2, 3,
        4, 5, 6,
        7, 8, 9,

        10, 11, 12,
        13, 14, 15,
        16, 17, 18,

        19, 20, 21,
        22, 23, 24,
        25, 26, 27
    ] );
    let f = e.multiply( &y );
    assert_eq!( f.rows(), 3 );
    assert_eq!( f.cols(), 3 );
    let g = &f * &y;
    assert_eq!( g.size(), 3 );
    let h = e.multiply_twice( &y );
    assert_eq!( h.size(), 3 );
    assert_eq!( g[0], h[0] );
    assert_eq!( g[1], h[1] );
    assert_eq!( g[2], h[2] );
}