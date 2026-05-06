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
fn full_non_consuming_addition() {
    let mut u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = &u + &v;
    assert_eq!( w[0], 7 );
    assert_eq!( w.size(), 10 );
    u[0] = 5;
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
fn full_non_consuming_subtraction() {
    let mut u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = &u - &v;
    assert_eq!( w[0], -1 );
    assert_eq!( w.size(), 10 );
    u[0] = 5;
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

#[test]
fn multiply_scalar() {
    let u = Vector::<f64>::new( 5, 3.0 );
    let v = 2.0 * u;
    assert_eq!( v[0], 6.0 );
}

#[test]
fn multiply_scalar_ref() {
    let u = Vector::<f64>::new( 5, 3.0 );
    let v = 2.0 * &u;
    assert_eq!( v[0], 6.0 );
}

#[test]
fn componentwise_multiply() {
    let u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    let w = &u * &v; // Full non-consuming componentwise multiplication
    for i in 0..w.size() {
        assert_eq!( w[i], 6.0 );
    }
    let x = u.clone() * &v; // Partial consuming componentwise multiplication
    for i in 0..x.size() {
        assert_eq!( x[i], 6.0 );
    }
    let y = u * v; // Consuming componentwise multiplication
    for i in 0..y.size() {
        assert_eq!( y[i], 6.0 );
    }
    let a = Vector::<i32>::create( vec![ 1, 2, 3, 4, 5 ] );
    let b = a.clone();
    let c = &a * &b;
    assert_eq!( c.size(), a.size() );
    for i in 0..c.size() {
        assert_eq!( c[i], ( i as i32 + 1 ) * ( i as i32 + 1 ) );
    }
}

#[test]
fn componentwise_divide() {
    let u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    let w = &u / &v; // Full non-consuming componentwise division
    for i in 0..w.size() {
        assert_eq!( w[i], 1.5 );
    }
    let x = u.clone() / &v; // Partial consuming componentwise division
    for i in 0..x.size() {
        assert_eq!( x[i], 1.5 );
    }
    let y = u / v; // Consuming componentwise division
    for i in 0..y.size() {
        assert_eq!( y[i], 1.5 );
    }
    let a = Vector::<i32>::create( vec![ 1, 2, 3, 4, 5 ] );
    let b = a.clone();
    let c = &a / &b;
    assert_eq!( c.size(), a.size() );
    for i in 0..c.size() {
        assert_eq!( c[i], 1 );
    }
}
