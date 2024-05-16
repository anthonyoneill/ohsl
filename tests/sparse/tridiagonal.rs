use ohsl::{Tridiagonal, Vector, Cmplx};

#[test]
fn unspecified_size_resize() {
    let mut mat = Tridiagonal::<f64>::empty();
    assert_eq!( mat.size(), 0 );
    assert_eq!( mat.subdiagonal().size(), 0 );
    assert_eq!( mat.maindiagonal().size(), 0 );
    assert_eq!( mat.superdiagonal().size(), 0 );
    mat.resize(5);
    assert_eq!( mat.size(), 5);
    assert_eq!( mat.subdiagonal().size(), 4 );
    assert_eq!( mat.maindiagonal().size(), 5 );
    assert_eq!( mat.superdiagonal().size(), 4 );
}

#[test]
fn specified_size(){
    let n = 5;
    let mat = Tridiagonal::<f64>::new(n);
    assert_eq!( mat.size(), n );
    assert_eq!( mat.subdiagonal().size(), n - 1 );
    assert_eq!( mat.maindiagonal().size(), n );
    assert_eq!( mat.superdiagonal().size(), n - 1 );
}

#[test]
fn constructor_with_vectors() {
    let n = 5;
    let sub = Vector::<f64>::new( n - 1, 1.0 );
    let main = Vector::<f64>::new( n, 2.0 );
    let sup = Vector::<f64>::new( n - 1, 3.0 );
    let mat = Tridiagonal::<f64>::with_vectors( sub, main, sup );
    assert_eq!( mat.size(), n );
    assert_eq!( mat.subdiagonal().size(), n - 1 );
    assert_eq!( mat.maindiagonal().size(), n );
    assert_eq!( mat.superdiagonal().size(), n - 1 );
}

#[test]
fn constructor_with_vecs() {
    let mat = Tridiagonal::<f64>::with_vecs( 
        vec![ 1.0, 2.0, 3.0 ], 
        vec![ 4.0, 5.0, 6.0, 7.0 ],
        vec![ 8.0, 9.0, 10.0 ] 
    );
    assert_eq!( mat.size(), 4 );
    assert_eq!( mat.subdiagonal().size(), 3 );
    assert_eq!( mat.maindiagonal().size(), 4 );
    assert_eq!( mat.superdiagonal().size(), 3 );
}

#[test]
fn constructor_with_elements() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 5);
    assert_eq!( mat.size(), 5);
    assert_eq!( mat.subdiagonal().size(), 4 );
    assert_eq!( mat.maindiagonal().size(), 5 );
    assert_eq!( mat.superdiagonal().size(), 4 );
}

#[test]
fn get_matrix_elements() {
    let mat = Tridiagonal::<f64>::with_vecs( 
        vec![ 1.0, 2.0 ], 
        vec![ 4.0, 5.0, 6.0 ],
        vec![ 8.0, 9.0 ] 
    );
    assert_eq!( mat[(0, 0)], 4.0 ); 
    assert_eq!( mat[(0, 1)], 8.0 );
    assert_eq!( mat[(1, 0)], 1.0 );
    assert_eq!( mat[(1, 1)], 5.0 );
    assert_eq!( mat[(1, 2)], 9.0 );
    assert_eq!( mat[(2, 1)], 2.0 );
    assert_eq!( mat[(2, 2)], 6.0 );
}

#[test]
#[should_panic]
fn not_an_element() {
    let mat = Tridiagonal::<f64>::with_vecs( 
        vec![ 1.0, 2.0 ], 
        vec![ 4.0, 5.0, 6.0 ],
        vec![ 8.0, 9.0 ] 
    );
    mat[(0, 2)];
}

#[test]
#[should_panic]
fn out_of_bounds() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    mat[(3, 0)];
}

#[test]
fn mutate_matrix_elements() {
    let mut mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    assert_eq!( mat[( 0, 0 )], 2.0 );
    mat[(0,0)] = 5.0;
    assert_eq!( mat[( 0, 0 )], 5.0 );
    mat[(0,0)] = 3.14;
    assert_eq!( mat[( 0, 0 )], 3.14 );
}

#[test]
fn clone_matrix() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = mat.clone();
    assert_eq!( mat2[( 0, 0 )], 2.0 );
    assert_eq!( mat2[( 0, 1 )], 3.0 );
    assert_eq!( mat2[( 1, 0 )], 1.0 );
    assert_eq!( mat2[( 1, 1 )], 2.0 );
    assert_eq!( mat2[( 1, 2 )], 3.0 );
    assert_eq!( mat2[( 2, 1 )], 1.0 );
    assert_eq!( mat2[( 2, 2 )], 2.0 );
}

#[test]
fn unary_negation() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = -mat;
    assert_eq!( mat2[( 0, 0 )], -2.0 );
    assert_eq!( mat2[( 0, 1 )], -3.0 );
    assert_eq!( mat2[( 1, 0 )], -1.0 );
    assert_eq!( mat2[( 1, 1 )], -2.0 );
    assert_eq!( mat2[( 1, 2 )], -3.0 );
    assert_eq!( mat2[( 2, 1 )], -1.0 );
    assert_eq!( mat2[( 2, 2 )], -2.0 );
}

#[test]
fn addition() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat3 = mat + mat2;
    assert_eq!( mat3[( 0, 0 )], 4.0 );
    assert_eq!( mat3[( 0, 1 )], 6.0 );
    assert_eq!( mat3[( 1, 0 )], 2.0 );
    assert_eq!( mat3[( 1, 1 )], 4.0 );
    assert_eq!( mat3[( 1, 2 )], 6.0 );
    assert_eq!( mat3[( 2, 1 )], 2.0 );
    assert_eq!( mat3[( 2, 2 )], 4.0 );
}

#[test]
fn subtraction() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat3 = mat - mat2;
    assert_eq!( mat3[( 0, 0 )], 0.0 );
    assert_eq!( mat3[( 0, 1 )], 0.0 );
    assert_eq!( mat3[( 1, 0 )], 0.0 );
    assert_eq!( mat3[( 1, 1 )], 0.0 );
    assert_eq!( mat3[( 1, 2 )], 0.0 );
    assert_eq!( mat3[( 2, 1 )], 0.0 );
    assert_eq!( mat3[( 2, 2 )], 0.0 );
}

#[test]
fn scalar_multiplication() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = mat.clone() * 2.0;
    assert_eq!( mat2[( 0, 0 )], 4.0 );
    assert_eq!( mat2[( 0, 1 )], 6.0 );
    assert_eq!( mat2[( 1, 0 )], 2.0 );
    assert_eq!( mat2[( 1, 1 )], 4.0 );
    assert_eq!( mat2[( 1, 2 )], 6.0 );
    assert_eq!( mat2[( 2, 1 )], 2.0 );
    assert_eq!( mat2[( 2, 2 )], 4.0 );
    let mat3 = 3.0 * mat;
    assert_eq!( mat3[( 0, 0 )], 6.0 );
    assert_eq!( mat3[( 0, 1 )], 9.0 );
    assert_eq!( mat3[( 1, 0 )], 3.0 );
    assert_eq!( mat3[( 1, 1 )], 6.0 );
    assert_eq!( mat3[( 1, 2 )], 9.0 );
    assert_eq!( mat3[( 2, 1 )], 3.0 );
    assert_eq!( mat3[( 2, 2 )], 6.0 );
}

#[test]
fn scalar_division() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let mat2 = mat / 2.0;
    assert_eq!( mat2[( 0, 0 )], 1.0 );
    assert_eq!( mat2[( 0, 1 )], 1.5 );
    assert_eq!( mat2[( 1, 0 )], 0.5 );
    assert_eq!( mat2[( 1, 1 )], 1.0 );
    assert_eq!( mat2[( 1, 2 )], 1.5 );
    assert_eq!( mat2[( 2, 1 )], 0.5 );
    assert_eq!( mat2[( 2, 2 )], 1.0 );
}

#[test]
fn add_assign() {
    let mut mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    mat += 2.0;
    assert_eq!( mat[( 0, 0 )], 4.0 );
    assert_eq!( mat[( 0, 1 )], 5.0 );
    assert_eq!( mat[( 1, 0 )], 3.0 );
    assert_eq!( mat[( 1, 1 )], 4.0 );
    assert_eq!( mat[( 1, 2 )], 5.0 );
    assert_eq!( mat[( 2, 1 )], 3.0 );
    assert_eq!( mat[( 2, 2 )], 4.0 );
}

#[test]
fn subtract_assign() {
    let mut mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    mat -= 2.0;
    assert_eq!( mat[( 0, 0 )], 0.0 );
    assert_eq!( mat[( 0, 1 )], 1.0 );
    assert_eq!( mat[( 1, 0 )], -1.0 );
    assert_eq!( mat[( 1, 1 )], 0.0 );
    assert_eq!( mat[( 1, 2 )], 1.0 );
    assert_eq!( mat[( 2, 1 )], -1.0 );
    assert_eq!( mat[( 2, 2 )], 0.0 );
}

#[test]
fn multiply_assign() {
    let mut mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    mat *= 3.0;
    assert_eq!( mat[( 0, 0 )], 6.0 );
    assert_eq!( mat[( 0, 1 )], 9.0 );
    assert_eq!( mat[( 1, 0 )], 3.0 );
    assert_eq!( mat[( 1, 1 )], 6.0 );
    assert_eq!( mat[( 1, 2 )], 9.0 );
    assert_eq!( mat[( 2, 1 )], 3.0 );
    assert_eq!( mat[( 2, 2 )], 6.0 );
}

#[test]
fn divide_assign() {
    let mut mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    mat /= 2.0;
    assert_eq!( mat[( 0, 0 )], 1.0 );
    assert_eq!( mat[( 0, 1 )], 1.5 );
    assert_eq!( mat[( 1, 0 )], 0.5 );
    assert_eq!( mat[( 1, 1 )], 1.0 );
    assert_eq!( mat[( 1, 2 )], 1.5 );
    assert_eq!( mat[( 2, 1 )], 0.5 );
    assert_eq!( mat[( 2, 2 )], 1.0 );
}

#[test]
fn vector_multiplication() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let vec = Vector::<f64>::ones( 3 );
    let ans = mat * vec;
    assert_eq!( ans[0], 5.0 );
    assert_eq!( ans[1], 6.0 );
    assert_eq!( ans[2], 3.0 );
}

#[test]
fn transpose() {
    let mut mat = Tridiagonal::<f64>::with_vecs( 
        vec![ 1.0, 2.0, 3.0 ], 
        vec![ 4.0, 5.0, 6.0, 7.0 ],
        vec![ 8.0, 9.0, 10.0 ] 
    );
    mat.transpose_in_place();
    assert_eq!( mat[( 0, 0 )], 4.0 );
    assert_eq!( mat[( 0, 1 )], 1.0 );
    assert_eq!( mat[( 1, 0 )], 8.0 );
    assert_eq!( mat[( 1, 1 )], 5.0 );
    assert_eq!( mat[( 1, 2 )], 2.0 );
    assert_eq!( mat[( 2, 1 )], 9.0 );
    assert_eq!( mat[( 2, 2 )], 6.0 );
    assert_eq!( mat[( 2, 3 )], 3.0 );
    assert_eq!( mat[( 3, 2 )], 10.0 );
    assert_eq!( mat[( 3, 3 )], 7.0 );
    let transpose = mat.transpose();
    assert_eq!( transpose[( 0, 0 )], 4.0 );
    assert_eq!( transpose[( 0, 1 )], 8.0 );
    assert_eq!( transpose[( 1, 0 )], 1.0 );
    assert_eq!( transpose[( 1, 1 )], 5.0 );
    assert_eq!( transpose[( 1, 2 )], 9.0 );
    assert_eq!( transpose[( 2, 1 )], 2.0 );
    assert_eq!( transpose[( 2, 2 )], 6.0 );
    assert_eq!( transpose[( 2, 3 )], 10.0 );
    assert_eq!( transpose[( 3, 2 )], 3.0 );
    assert_eq!( transpose[( 3, 3 )], 7.0 );
}

#[test]
fn determinant() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    assert_eq!( mat.det(), -4.0 );
}

#[test]
fn convert_to_dense() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let dense = mat.convert();
    assert_eq!( dense[(0,0)], 2.0 );
    assert_eq!( dense[(0,1)], 3.0 );
    assert_eq!( dense[(0,2)], 0.0 );
    assert_eq!( dense[(1,0)], 1.0 );
    assert_eq!( dense[(1,1)], 2.0 );
    assert_eq!( dense[(1,2)], 3.0 );
    assert_eq!( dense[(2,0)], 0.0 );
    assert_eq!( dense[(2,1)], 1.0 );
    assert_eq!( dense[(2,2)], 2.0 );
    let mat2 = Tridiagonal::<f64>::with_vecs( 
        vec![ 1.0 ], 
        vec![ 2.0, 3.0 ],
        vec![ 4.0 ] 
    );
    let dense2 = mat2.convert();
    assert_eq!( dense2[(0,0)], 2.0 );
    assert_eq!( dense2[(0,1)], 4.0 );
    assert_eq!( dense2[(1,0)], 1.0 );
    assert_eq!( dense2[(1,1)], 3.0 );
}

#[test]
fn conjugate() {
    let mat = Tridiagonal::<Cmplx>::with_elements( 
        Cmplx::new( 1.0, 1.0 ), 
        Cmplx::new( 2.0, 2.0 ), 
        Cmplx::new( 3.0, 3.0 ), 
    3);
    let conj = mat.conj();
    assert_eq!( conj[( 0, 0 )], Cmplx::new( 2.0, -2.0 ) );
    assert_eq!( conj[( 0, 1 )], Cmplx::new( 3.0, -3.0 ) );
    assert_eq!( conj[( 1, 0 )], Cmplx::new( 1.0, -1.0 ) );
    assert_eq!( conj[( 1, 1 )], Cmplx::new( 2.0, -2.0 ) );
    assert_eq!( conj[( 1, 2 )], Cmplx::new( 3.0, -3.0 ) );
    assert_eq!( conj[( 2, 1 )], Cmplx::new( 1.0, -1.0 ) );
    assert_eq!( conj[( 2, 2 )], Cmplx::new( 2.0, -2.0 ) );
}

#[test]
fn solve() {
    let mat = Tridiagonal::<f64>::with_elements( 1.0, 2.0, 3.0, 3);
    let b = Vector::<f64>::ones( 3 );
    let x = mat.solve( &b );
    assert_eq!( x[0], -1.0 );
    assert_eq!( x[1], 1.0 );
    assert_eq!( x[2], 0.0 );
} 