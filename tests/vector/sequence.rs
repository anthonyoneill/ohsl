use ohsl::vector::Vector;

#[test]
fn linspace() {
    let v = Vector::<f64>::linspace( 0.0, 1.0, 11 );
    assert_eq!( v[10], 1.0 );
    assert_eq!( v[0], 0.0 );
    assert_eq!( v[5], 0.5 );
}

#[test]
fn powspace() {
    let v = Vector::<f64>::powspace( 0.0, 1.0, 11, 2.0 );
    assert_eq!( v[10], 1.0 );
    assert_eq!( v[0], 0.0 );
    assert_eq!( v[5], 0.25 );
}