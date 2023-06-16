use ohsl::vector::{Vector, Vec64};

#[test]
fn unspecified_size() {
    let v = Vector::<i32>::empty();
    assert_eq!( v.size(), 0 );
}

#[test]
fn specified_size() {
    let mut v = Vec64::new( 10, 3.14 );
    assert_eq!( v.size(), 10 );
    assert_eq!( v[3], 3.14 );
    v[3] = 7.0;
    assert_eq!( v[3], 7.0 );
}

#[test]
fn from_vec() {
    let v = Vector::<u32>::create( vec![ 1, 2, 3 ] );
    assert_eq!( v[1], 2 );
    assert_eq!( v.size(), 3 );
    let w = v.clone();
    assert_eq!( w[2], 3 );
    assert_eq!( w.size(), 3 );
}

#[test]
fn assigment() {
    let w = Vector::<i32>::new( 10, 4 );
    let v = w;
    assert_eq!( v[0], 4 );
    assert_eq!( v.size(), 10 );
}

#[test]
fn zeros() {
    let zeros = Vector::<i32>::zeros( 5 );
    for i in 0..5 {
        assert_eq!( zeros[i], 0 );
    }
}

#[test]
fn ones() {
    let ones = Vector::<f64>::ones( 7 );
    for i in 0..7 {
        assert_eq!( ones[i], 1. );
    }
}