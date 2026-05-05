use ohsl::triad::Triad;

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