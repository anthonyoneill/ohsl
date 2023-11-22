use ohsl::Polynomial;

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
}