use ohsl::{Banded, Vec64, Cmplx, Zero, Vector};

#[test]
fn empty_constructor() {
    let mut mat = Banded::<f64>::empty();
    assert_eq!( mat.size(), 0 );
    assert_eq!( mat.size_below(), 0 );
    assert_eq!( mat.size_above(), 0 );
    assert_eq!( mat.compact().rows(), 0 );
    assert_eq!( mat.compact().cols(), 0 );
    mat.resize( 3, 1, 1 );
    assert_eq!( mat.size(), 3 );
    assert_eq!( mat.size_below(), 1 );
    assert_eq!( mat.size_above(), 1 );
    assert_eq!( mat.compact().rows(), 3 );
    assert_eq!( mat.compact().cols(), 1 + 1 + 1 );
}

#[test]
fn constructor() {
    let mat = Banded::<f64>::new( 5, 2, 1, 0.0 );
    assert_eq!( mat.size(), 5 );
    assert_eq!( mat.size_below(), 2 );
    assert_eq!( mat.size_above(), 1 );
    assert_eq!( mat.compact().rows(), 5 );
    assert_eq!( mat.compact().cols(), 2 + 1 + 1 );
}

#[test]
fn indexing() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    assert_eq!( mat[(0,0)], 1.0 );
    assert_eq!( mat[(0,1)], 1.0 );
    assert_eq!( mat[(1,0)], 1.0 );
    assert_eq!( mat[(1,1)], 1.0 );
    assert_eq!( mat[(1,2)], 1.0 );
    assert_eq!( mat[(2,1)], 1.0 );
    assert_eq!( mat[(2,2)], 1.0 );
    mat[(0,0)] = 2.0;
    mat[(0,1)] = 3.0;
    assert_eq!( mat[(0,0)], 2.0 );
    assert_eq!( mat[(0,1)], 3.0 );
}

#[test]
#[should_panic]
fn indexing_out_of_bounds() {
    let mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    let _ = mat[(3,0)];
}

#[test]
#[should_panic]
fn indexing_out_of_band() {
    let mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    let _ = mat[(0,2)];
}

#[test]
fn fill_band() {
    let mut mat = Banded::<f64>::new( 5, 2, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    for i in 0..5 { assert_eq!( mat[(i,i)], 3.0 ); }
    for i in 0..4 { assert_eq!( mat[(i,i+1)], 4.0 ); }
    for i in 1..5 { assert_eq!( mat[(i,i-1)], 2.0 ); }
    for i in 2..4 { assert_eq!( mat[(i,i-2)], 1.0 ); }
}

#[test]
#[should_panic]
fn fill_band_out_of_bounds() {
    let mut mat = Banded::<f64>::new( 5, 2, 1, 1.0 );
    mat.fill_band( -3, 2.0 );
}

#[test]
fn display() {
    let mut mat = Banded::<f64>::new( 5, 2, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let display = format!("{}", mat);
    let mut string = String::new();
    string.push_str( "\t3\t4\t*\t*\t*\n" );
    string.push_str( "\t2\t3\t4\t*\t*\n" );
    string.push_str( "\t1\t2\t3\t4\t*\n" );
    string.push_str( "\t*\t1\t2\t3\t4\n" );
    string.push_str( "\t*\t*\t1\t2\t3\n" );
    string.push_str( "\n" );
    assert_eq!( display, string );
}

#[test]
fn clone() {
    let mut mat = Banded::<f64>::new( 5, 2, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = mat.clone();
    assert_eq!( mat, mat2 );
}

#[test]
fn negation() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    println!("mat = {}", mat);
    let mat2 = -&mat; // Non-consuming negation
    println!("mat2 = {}", mat2);
    assert_eq!( mat2[(0,0)], -3.0 );
    assert_eq!( mat2[(0,1)], -4.0 );
    assert_eq!( mat2[(1,0)], -2.0 );
    assert_eq!( mat2[(1,1)], -3.0 );
    assert_eq!( mat2[(1,2)], -4.0 );
    assert_eq!( mat2[(2,1)], -2.0 );
    assert_eq!( mat2[(2,2)], -3.0 );
    let mat3 = -mat2; // Consuming negation
    assert_eq!( mat3[(0,0)], 3.0 );
    assert_eq!( mat3[(0,1)], 4.0 );
    assert_eq!( mat3[(1,0)], 2.0 );
    assert_eq!( mat3[(1,1)], 3.0 );
    assert_eq!( mat3[(1,2)], 4.0 );
    assert_eq!( mat3[(2,1)], 2.0 );
    assert_eq!( mat3[(2,2)], 3.0 );
}

#[test]
fn addition() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = mat.clone();
    let mat3 = mat + mat2; // Consuming addition
    assert_eq!( mat3[(0,0)], 6.0 );
    assert_eq!( mat3[(0,1)], 8.0 );
    assert_eq!( mat3[(1,0)], 4.0 );
    assert_eq!( mat3[(1,1)], 6.0 );
    assert_eq!( mat3[(1,2)], 8.0 );
    assert_eq!( mat3[(2,1)], 4.0 );
    assert_eq!( mat3[(2,2)], 6.0 );
    let mat4 = &mat3 + &mat3; // Non-consuming addition
    assert_eq!( mat4[(0,0)], 12.0 );
    assert_eq!( mat4[(0,1)], 16.0 );
    assert_eq!( mat4[(1,0)], 8.0 );
    assert_eq!( mat4[(1,1)], 12.0 );
    assert_eq!( mat4[(1,2)], 16.0 );
    assert_eq!( mat4[(2,1)], 8.0 );
    assert_eq!( mat4[(2,2)], 12.0 );
}

#[test]
fn subtraction() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = mat.clone();
    let mat3 = &mat - &mat2; // Non-Consuming subtraction
    assert_eq!( mat3[(0,0)], 0.0 );
    assert_eq!( mat3[(0,1)], 0.0 );
    assert_eq!( mat3[(1,0)], 0.0 );
    assert_eq!( mat3[(1,1)], 0.0 );
    assert_eq!( mat3[(1,2)], 0.0 );
    assert_eq!( mat3[(2,1)], 0.0 );
    assert_eq!( mat3[(2,2)], 0.0 );
    let mat4 = mat3 - mat; // Consuming subtraction
    assert_eq!( mat4[(0,0)], -3.0 );
    assert_eq!( mat4[(0,1)], -4.0 );
    assert_eq!( mat4[(1,0)], -2.0 );
    assert_eq!( mat4[(1,1)], -3.0 );
    assert_eq!( mat4[(1,2)], -4.0 );
    assert_eq!( mat4[(2,1)], -2.0 );
    assert_eq!( mat4[(2,2)], -3.0 );
}

#[test]
fn scalar_multiplication() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = &mat * 2.0; // Non-consuming scalar multiplication
    assert_eq!( mat2[(0,0)], 6.0 );
    assert_eq!( mat2[(0,1)], 8.0 );
    assert_eq!( mat2[(1,0)], 4.0 );
    assert_eq!( mat2[(1,1)], 6.0 );
    assert_eq!( mat2[(1,2)], 8.0 );
    assert_eq!( mat2[(2,1)], 4.0 );
    assert_eq!( mat2[(2,2)], 6.0 );
    let mat3 = mat * 4.0; // Consuming scalar multiplication
    assert_eq!( mat3[(0,0)], 12.0 );
    assert_eq!( mat3[(0,1)], 16.0 );
    assert_eq!( mat3[(1,0)], 8.0 );
    assert_eq!( mat3[(1,1)], 12.0 );
    assert_eq!( mat3[(1,2)], 16.0 );
    assert_eq!( mat3[(2,1)], 8.0 );
    assert_eq!( mat3[(2,2)], 12.0 );
}

#[test]
fn scalar_division() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = &mat / 2.0; // Non-consuming scalar division
    assert_eq!( mat2[(0,0)], 1.5 );
    assert_eq!( mat2[(0,1)], 2.0 );
    assert_eq!( mat2[(1,0)], 1.0 );
    assert_eq!( mat2[(1,1)], 1.5 );
    assert_eq!( mat2[(1,2)], 2.0 );
    assert_eq!( mat2[(2,1)], 1.0 );
    assert_eq!( mat2[(2,2)], 1.5 );
    let mat3 = mat / 4.0; // Consuming scalar division
    assert_eq!( mat3[(0,0)], 0.75 );
    assert_eq!( mat3[(0,1)], 1.0 );
    assert_eq!( mat3[(1,0)], 0.5 );
    assert_eq!( mat3[(1,1)], 0.75 );
    assert_eq!( mat3[(1,2)], 1.0 );
    assert_eq!( mat3[(2,1)], 0.5 );
    assert_eq!( mat3[(2,2)], 0.75 );
}

#[test]
fn addition_assignment() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = mat.clone();
    mat += &mat2; // Non-consuming addition assignment
    assert_eq!( mat[(0,0)], 6.0 );
    assert_eq!( mat[(0,1)], 8.0 );
    assert_eq!( mat[(1,0)], 4.0 );
    assert_eq!( mat[(1,1)], 6.0 );
    assert_eq!( mat[(1,2)], 8.0 );
    assert_eq!( mat[(2,1)], 4.0 );
    assert_eq!( mat[(2,2)], 6.0 );
    mat += mat2; // Consuming addition assignment
    assert_eq!( mat[(0,0)], 9.0 );
    assert_eq!( mat[(0,1)], 12.0 );
    assert_eq!( mat[(1,0)], 6.0 );
    assert_eq!( mat[(1,1)], 9.0 );
    assert_eq!( mat[(1,2)], 12.0 );
    assert_eq!( mat[(2,1)], 6.0 );
    assert_eq!( mat[(2,2)], 9.0 );
}

#[test]
fn subtraction_assignment() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let mat2 = mat.clone();
    mat -= &mat2; // Non-consuming subtraction assignment
    assert_eq!( mat[(0,0)], 0.0 );
    assert_eq!( mat[(0,1)], 0.0 );
    assert_eq!( mat[(1,0)], 0.0 );
    assert_eq!( mat[(1,1)], 0.0 );
    assert_eq!( mat[(1,2)], 0.0 );
    assert_eq!( mat[(2,1)], 0.0 );
    assert_eq!( mat[(2,2)], 0.0 );
    mat -= mat2; // Consuming subtraction assignment
    assert_eq!( mat[(0,0)], -3.0 );
    assert_eq!( mat[(0,1)], -4.0 );
    assert_eq!( mat[(1,0)], -2.0 );
    assert_eq!( mat[(1,1)], -3.0 );
    assert_eq!( mat[(1,2)], -4.0 );
    assert_eq!( mat[(2,1)], -2.0 );
    assert_eq!( mat[(2,2)], -3.0 );
}

#[test]
fn scalar_multiplication_assignement() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 1.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    mat *= 2.0; 
    assert_eq!( mat[(0,0)], 6.0 );
    assert_eq!( mat[(0,1)], 8.0 );
    assert_eq!( mat[(1,0)], 4.0 );
    assert_eq!( mat[(1,1)], 6.0 );
    assert_eq!( mat[(1,2)], 8.0 );
    assert_eq!( mat[(2,1)], 4.0 );
    assert_eq!( mat[(2,2)], 6.0 );
}

#[test]
fn scalar_division_assignment() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    mat /= 2.0; 
    assert_eq!( mat[(0,0)], 1.5 );
    assert_eq!( mat[(0,1)], 2.0 );
    assert_eq!( mat[(1,0)], 1.0 );
    assert_eq!( mat[(1,1)], 1.5 );
    assert_eq!( mat[(1,2)], 2.0 );
    assert_eq!( mat[(2,1)], 1.0 );
    assert_eq!( mat[(2,2)], 1.5 );
}

#[test]
fn constant_addition_assignment() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    mat += 2.0; 
    assert_eq!( mat[(0,0)], 5.0 );
    assert_eq!( mat[(0,1)], 6.0 );
    assert_eq!( mat[(1,0)], 4.0 );
    assert_eq!( mat[(1,1)], 5.0 );
    assert_eq!( mat[(1,2)], 6.0 );
    assert_eq!( mat[(2,1)], 4.0 );
    assert_eq!( mat[(2,2)], 5.0 );
}

#[test]
fn constant_subtraction_assignment() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    mat -= 2.0; 
    assert_eq!( mat[(0,0)], 1.0 );
    assert_eq!( mat[(0,1)], 2.0 );
    assert_eq!( mat[(1,0)], 0.0 );
    assert_eq!( mat[(1,1)], 1.0 );
    assert_eq!( mat[(1,2)], 2.0 );
    assert_eq!( mat[(2,1)], 0.0 );
    assert_eq!( mat[(2,2)], 1.0 );
}

#[test]
fn matrix_vector_multiplication() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    println!("mat = {}", mat);
    let v = Vec64::new( 3, 1.0 );
    println!("v = {}", v);
    let w = &mat * &v; // Non-consuming matrix-vector multiplication
    assert_eq!( w[0], 7.0 );
    assert_eq!( w[1], 9.0 );
    assert_eq!( w[2], 5.0 );
    let u = Vec64::create( vec![1.0, -1.0, 2.0] );
    let w2 = mat * u; // Consuming matrix-vector multiplication
    assert_eq!( w2[0], -1.0 );
    assert_eq!( w2[1],  7.0 );
    assert_eq!( w2[2],  4.0 );
}

#[test]
fn determinant() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    assert_eq!( mat.det(), -21.0 );
}

#[test]
fn complex_determinant() {
    let mut mat = Banded::<Cmplx>::new( 3, 1, 1, Cmplx::zero() );
    mat.fill_band( 0, Cmplx::new( 3.0, 1.0 ) );
    mat.fill_band( -1, Cmplx::new( 2.0, 0.0 ) );
    mat.fill_band( 1, Cmplx::new( 4.0, 0.0 ) );
    assert_eq!( mat.det(), Cmplx::new( -30.0, 10.0 ) );
}

#[test]
fn solve() {
    let mut mat = Banded::<f64>::new( 3, 1, 1, 2.0 );
    mat.fill_band( 0, 3.0 );
    mat.fill_band( -1, 2.0 );
    mat.fill_band( 1, 4.0 );
    let v = Vec64::new( 3, 1.0 );
    let w = mat.solve( &v );
    assert!( (w[0] - (-5.0 / 21.0)).abs() < f64::EPSILON ); 
    assert!( (w[1] - ( 9.0 / 21.0)).abs() < f64::EPSILON );
    assert!( (w[2] - ( 1.0 / 21.0)).abs() < f64::EPSILON );
}

#[test]
fn complex_solve() {
    let mut mat = Banded::<Cmplx>::new( 3, 1, 1, Cmplx::zero() );
    mat.fill_band( 0, Cmplx::new( 3.0, 1.0 ) );
    mat.fill_band( -1, Cmplx::new( 2.0, 0.0 ) );
    mat.fill_band( 1, Cmplx::new( 4.0, 0.0 ) );
    let v = Vector::<Cmplx>::new( 3, Cmplx::new( 1.0, 0.0 ) );
    let w = mat.solve( &v );
    let factor = Cmplx::new( 0.1, 0.1 );
    let w0 = factor * Cmplx::new( -1.0, 0.0 );
    let w1 = factor * Cmplx::new(  2.0, -1.0 );
    let w2 = factor * Cmplx::new(  0.0, -1.0 );
    assert!( (w[0] - w0).abs() < f64::EPSILON ); 
    assert!( (w[1] - w1).abs() < f64::EPSILON );
    assert!( (w[2] - w2).abs() < f64::EPSILON );
}