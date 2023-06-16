use ohsl::vector::Vector;

#[test]
fn norm_1() {
    let u = Vector::<f64>::create( vec![ 1.0, -2.0, 2.0, 3.0, 5.0 ] );
    let l1 = u.norm_1();
    assert_eq!( l1, 13.0 );
}

#[test]
fn norm_2() {
    let v = Vector::<f64>::create( vec![ 3.0, 4.0 ] );
    let l2 = v.norm_2();
    assert_eq!( l2, 5.0 );
}

#[test]
fn norm_p() {
    let v = Vector::<f64>::create( vec![ 3.0, 4.0 ] );
    let lp = v.norm_p( 2.0 );
    assert_eq!( lp, 5.0 );
    let lp = v.norm_p( 3.0 );
    assert_eq!( lp, 4.497941445275415 );
    let lp = v.norm_p( 2.5 );
    assert_eq!( lp, 4.688140842343588 );
}

#[test]
fn norm_inf() {
    let u = Vector::<f64>::create( vec![ 1.0, -2.0, 2.0, 3.0, 5.0 ] );
    let linf = u.norm_inf();
    assert_eq!( linf, 5.0 );
}