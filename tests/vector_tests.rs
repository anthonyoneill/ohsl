use ohsl::vector::{Vector, Vec64};
use ohsl::complex::Cmplx;

#[test]
fn test_vector_unspecified_size() {
    let v = Vector::<i32>::empty();
    assert_eq!( v.size(), 0 );
}

#[test]
fn test_vector_specified_size() {
    let mut v = Vec64::new( 10, 3.14 );
    assert_eq!( v.size(), 10 );
    assert_eq!( v[3], 3.14 );
    v[3] = 7.0;
    assert_eq!( v[3], 7.0 );
}

#[test]
fn test_initialise_from_vec() {
    let v = Vector::<u32>::create( vec![ 1, 2, 3 ] );
    assert_eq!( v[1], 2 );
    assert_eq!( v.size(), 3 );
    let w = v.clone();
    assert_eq!( w[2], 3 );
    assert_eq!( w.size(), 3 );
}

#[test]
fn test_vector_assigment() {
    let w = Vector::<i32>::new( 10, 4 );
    let v = w;
    assert_eq!( v[0], 4 );
    assert_eq!( v.size(), 10 );
}

#[test]
fn test_vector_unary_minus() {
    let v = Vector::<i32>::new( 10, 4 );
    let w = -v;
    assert_eq!( w[0], -4 ); 
}

#[test]
fn test_vector_binary_addition() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u + v;
    assert_eq!( w[0], 7 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn test_vector_binary_subtraction() {
    let u = Vector::<i32>::new( 10, 3 );
    let v = Vector::<i32>::new( 10, 4 );
    let w = u - v;
    assert_eq!( w[0], -1 );
    assert_eq!( w.size(), 10 );
}

#[test]
fn test_vector_addassign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    u += v;
    assert_eq!( u[0], 5.0 );
    u += 5.0;
    assert_eq!( u[1], 10.0 );
}

#[test]
fn test_vector_subassign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    let v = Vector::<f64>::new( 5, 2.0 );
    u -= v;
    assert_eq!( u[0], 1.0 );
    assert_eq!( u.size(), 5 );
    u -= 1.0;
    assert_eq!( u[4], 0.0 );
}

#[test]
fn test_vector_mulassign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u *= 2.0;
    assert_eq!( u[0], 6.0 );
}

#[test]
fn test_vector_divassign() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u /= 2.0;
    assert_eq!( u[0], 1.5 );
}

#[test]
fn test_resize_assign_clear() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.resize( 10 );
    assert_eq!( u[0], 3.0 );
    assert_eq!( u[9], 0.0 );
    u.assign( 7.0 );
    assert_eq!( u[3], 7.0 );
    u.clear();
    assert_eq!( u.size(), 0 );
}

#[test]
fn test_swap() {
    let mut u = Vector::<f64>::new( 4, 1.0 );
    u[0] = 7.0;
    u.swap( 0, 2 );
    assert_eq!( u[2], 7.0 );
}

#[test]
fn test_dot_product() {
    let u = Vector::<u32>::create( vec![ 1, 2, 3 ] );
    let v = Vector::<u32>::create( vec![ 4, 5, 6 ] );
    let dot = u.dot( v );
    assert_eq!( dot, 32 );
}

#[test]
fn test_sum_product() {
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

/*#[test]
fn test_vector_file_output() {
    let u = Vector::<u32>::create( vec![ 1, 2, 3, 4, 5 ] );
    u.output( "./output.txt" );
}*/

#[test]
fn test_push_pop() {
    let mut u = Vector::<f64>::new( 5, 3.0 );
    u.push( 7.0 );
    assert_eq!( u[5], 7.0 );
    assert_eq!( u.size(), 6 );
    let v = u.pop();
    assert_eq!( v, 7.0 );
    assert_eq!( u.size(), 5 );
}

#[test]
fn test_vector_find() {
    let u = Vector::<u32>::create( vec![ 1, 2, 7, 4, 5 ] );
    let index = u.find( 7 );
    assert_eq!( index, 2 );
    //let index = u.find( 8 );
}

#[test]
fn test_vector_insert() {
    let mut u = Vector::<u32>::create( vec![ 1, 2, 7, 4, 5 ] );
    u.insert( 0, 8 );
    assert_eq!( u[0], 8 );
    assert_eq!( u[1], 1 );
}

#[test]
fn test_vector_abs() {
    let u = Vector::<i32>::new( 5, -2 );
    let v = u.abs();
    assert_eq!( v[2], 2 );
}

#[test]
fn test_vector_conjugate() {
    let v = Vector::<Cmplx>::new( 5, Cmplx::new( 1.0, 2.0 ) );
    assert_eq!( v[0].real, 1.0 );
    let conj = v.conj();
    assert_eq!( conj[0].imag, -2.0 );
}

#[test]
fn test_vector_real() {
    let v = Vector::<Cmplx>::new( 5, Cmplx::new( 1.0, 2.0 ) );
    let real = v.real();
    assert_eq!( real[0], 1.0 );
}

#[test]
fn test_linspace() {
    let v = Vector::<f64>::linspace( 0.0, 1.0, 11 );
    assert_eq!( v[10], 1.0 );
    assert_eq!( v[0], 0.0 );
    assert_eq!( v[5], 0.5 );
}

#[test]
fn test_powspace() {
    let v = Vector::<f64>::powspace( 0.0, 1.0, 11, 2.0 );
    assert_eq!( v[10], 1.0 );
    assert_eq!( v[0], 0.0 );
    assert_eq!( v[5], 0.25 );
}

#[test]
fn test_norms() {
    let u = Vector::<f64>::create( vec![ 1.0, -2.0, 2.0, 3.0, 5.0 ] );
    let l1 = u.norm_1();
    assert_eq!( l1, 13.0 );
    let v = Vector::<f64>::create( vec![ 3.0, 4.0 ] );
    let l2 = v.norm_2();
    assert_eq!( l2, 5.0 );
    let lp = v.norm_p( 2.0 );
    assert_eq!( lp, 5.0 );
    let linf = u.norm_inf();
    assert_eq!( linf, 5.0 );
}

#[test]
fn test_vector_random() {
    let v = Vec64::random( 5 );
    assert_eq!( v.size(), 5 );
    assert!( v[0] > 0.0 );
    assert!( v[0] < 1.0 );
    assert!( v[0] - v[1] != 0.0 );
}

#[test]
fn test_zeros_ones() {
    let zeros = Vector::<i32>::zeros( 5 );
    assert_eq!( zeros[4], 0 );
    let ones = Vector::<f64>::ones( 7 );
    assert_eq!( ones[6], 1.0 );
}

#[test]
fn test_sort() {
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