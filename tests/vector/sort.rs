use ohsl::vector::{Vector, Vec64};

#[test]
fn sort() {
    let mut u = Vector::<u32>::create( vec![ 5, 4, 3, 2, 1 ] );
    u.sort();
    assert_eq!( u[0], 1 );
    assert_eq!( u[1], 2 );
    assert_eq!( u[2], 3 );
    assert_eq!( u[3], 4 );
    assert_eq!( u[4], 5 );
    let mut v = Vec64::create( vec![ 5.0, 4.0, 3.0, 2.0, 1.0 ] );
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!( v[0], 1.0 );
    assert_eq!( v[1], 2.0 );
    assert_eq!( v[2], 3.0 );
    assert_eq!( v[3], 4.0 );
    assert_eq!( v[4], 5.0 );
}