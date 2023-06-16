use ohsl::vector::Vec64;

#[test]
fn test_vector_random() {
    let v = Vec64::random( 5 );
    assert_eq!( v.size(), 5 );
    assert!( v[0] > 0.0 );
    assert!( v[0] < 1.0 );
    assert!( v[0] - v[1] != 0.0 );
}