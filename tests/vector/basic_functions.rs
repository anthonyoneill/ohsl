use ohsl::vector::{Vector, Vec64};

#[test]
fn resize() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.resize( 10 );
    assert_eq!( u[0], 3.0 );
    assert_eq!( u[9], 0.0 );
}

#[test]
fn assign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.assign( 7.0 );
    assert_eq!( u[3], 7.0 );
}

#[test]
fn clear() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.clear();
    assert_eq!( u.size(), 0 );
}

#[test]
fn swap() {
    let mut u = Vector::<f64>::new( 4, 1.0 );
    u[0] = 7.0;
    u.swap( 0, 2 );
    assert_eq!( u[2], 7.0 );
}

#[test]
fn dot_product() {
    let u = Vector::<u32>::create( vec![ 1, 2, 3 ] );
    let v = Vector::<u32>::create( vec![ 4, 5, 6 ] );
    let dot = u.dot( &v );
    assert_eq!( dot, 32 );
    let u = Vec64::create( vec![ 1.0, 2.0, 3.0 ] );
    let v = Vec64::create( vec![ 4.0, 5.0, 6.0 ] );
    let dotf64 = u.dot( &v );
    assert_eq!( dotf64, 32.0 );
}

#[test]
fn sum_product() {
    let u = Vector::<u32>::create( vec![ 1, 2, 3, 4, 5 ] );
    let sum = u.sum();
    assert_eq!( sum, 15 );
    let sum_slice = u.sum_slice( 1, 3 );
    assert_eq!( sum_slice, 9 );
    let product = u.product();
    assert_eq!( product, 120 );
    let product_slice = u.product_slice( 2, 3 );
    assert_eq!( product_slice, 12 );
}

#[test]
fn push() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.push( 7.0 );
    assert_eq!( u[5], 7.0 );
    assert_eq!( u.size(), 6 );
}

#[test]
fn pop() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    let v = u.pop();
    assert_eq!( v, 3.0 );
    assert_eq!( u.size(), 4 );
}

#[test]
fn find() {
    let u = Vector::<u32>::create( vec![ 1, 2, 7, 4, 5 ] );
    let index = u.find( 7 );
    assert_eq!( index, 2 );
}

#[test]
fn insert() {
    let mut u = Vector::<u32>::create( vec![ 1, 2, 7, 4, 5 ] );
    u.insert( 0, 8 );
    assert_eq!( u[0], 8 );
    assert_eq!( u[1], 1 );
}

#[test]
fn abs() {
    let u = Vector::<i32>::new( 5, -2 );
    let v = u.abs();
    assert_eq!( v[2], 2 );
}