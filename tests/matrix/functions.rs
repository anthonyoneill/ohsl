use ohsl::matrix::Mat64;

#[test]
fn norms() {
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

#[test]
fn determinant() {
    let mut a = Mat64::new( 2, 2, 4.0 );
    a[0][1] = 3.0;
    a[1][0] = 6.0; 
    a[1][1] = 3.0;
    let det = a.determinant();
    assert_eq!( det, -6.0 );
}

#[test]
fn inverse() {
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