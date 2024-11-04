use ohsl::vector::Vector;

#[test]
fn unary_minus() {
    let v = Vector::<i32>::new( 10, 4 );
    let w = -v;
    assert_eq!( w[0], -4 ); 
}

#[test]
fn binary_addition() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u + v;
    assert_eq!( w[0], 7 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn non_consuming_addition() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u + &v;
    assert_eq!( w[0], 7 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn binary_subtraction() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u - v;
    assert_eq!( w[0], -1 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn non_consuming_subtraction() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u - &v;
    assert_eq!( w[0], -1 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn add_assign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    u += v;
    assert_eq!( u[0], 5.0 );
    u += 5.0;
    assert_eq!( u[1], 10.0 );
}

#[test]
fn subtract_assign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    u -= v;
    assert_eq!( u[0], 1.0 );
    assert_eq!( u.size(), 5 );
    u -= 1.0;
    assert_eq!( u[4], 0.0 );
}

#[test]
fn multiply_assign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u *= 2.0;
    assert_eq!( u[0], 6.0 );
}

#[test]
fn divide_assign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u /= 2.0;
    assert_eq!( u[0], 1.5 );
}