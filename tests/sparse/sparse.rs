mod tridiagonal;
mod banded;

use ohsl::sparse_matrix::{Sparse, MatrixState, Axis};
use ohsl::vector::{Vector, Vec64};

#[test]
fn create_matrix() {
    let m = Sparse::new();
    assert_eq!(m.state(), MatrixState::CREATED);
    assert_eq!(m.diag(), vec![]);
}

#[test]
fn add_element() {
    let mut s = Sparse::new();

    s.add_element(0, 0, 1.0);
    assert_eq!(s.num_rows(), 1);
    assert_eq!(s.num_cols(), 1);
    assert_eq!(s.size(), (1, 1));
    assert_eq!(s.diag().len(), 1);

    s.add_element(100, 100, 1.0);
    assert_eq!(s.num_rows(), 101);
    assert_eq!(s.num_cols(), 101);
    assert_eq!(s.size(), (101, 101));
    assert_eq!(s.diag().len(), 101);
}

#[test]
fn get() {
    let mut s = Sparse::new();
    s.add_element(0, 0, 1.0);
    assert_eq!( s.get(0, 0).unwrap(), 1.0 );
}

#[test]
fn identity() {
    // Check identity matrices of each (small) size
    for k in 1..10 {
        let ik = Sparse::identity(k);

        // Basic size checks
        assert_eq!(ik.num_rows(), k);
        assert_eq!(ik.num_cols(), k);
        assert_eq!(ik.size(), (k, k));
        assert_eq!(ik.elements().len(), k);

        for v in 0..k {
            // Check each row/ col head is the same element, and this element is on the diagonal
            let ro = ik.hdr(Axis::ROWS, v).unwrap();
            let co = ik.hdr(Axis::COLS, v).unwrap();
            let d0 = ik.get(v, v).unwrap();
            assert_eq!(ro, co);
            assert_eq!(ik[ro].val(), d0);
            assert_eq!(ik[co].val(), d0);
        }
    }
}

#[test]
fn solve() {
    let mut s = Sparse::from_triplets(vec![
        (0, 0, 1.0),
        (0, 1, 2.0),
        (1, 0, 3.0),
        (1, 1, 4.0),
        (1, 2, 5.0),
        (2, 1, 6.0),
        (2, 2, 7.0),
    ]);

    // ohsl Vector solve
    let b = Vec64::new( 3, 1.0 );
    let solution = s.solve( b ).unwrap();
    assert_eq!( solution[0], 6.0/44.0 );
    assert_eq!( solution[1], 19.0/44.0 );
    assert_eq!( solution[2], -10.0/44.0 );
    
    // std Vec _solve
    let rhs = vec![1.0, 1.0, 1.0];
    let soln = s._solve(rhs).unwrap();
    let correct = vec![6.0/44.0, 19.0/44.0, -10.0/44.0];
    assert_eq!( soln, correct );
}

/* Sparse2 */

use ohsl::sparse2::Sparse2;

fn test_sparse_matrix() -> Sparse2::<f64> {
    let val = vec![ 3.0, 4.0, 7.0, 1.0, 5.0, 2.0, 9.0, 6.0, 5.0 ];
    let row_index = vec![ 0, 1, 2, 0, 2, 0, 2, 4, 4 ];
    let col_start = vec![ 0, 1, 3, 5, 8, 9 ];
    Sparse2::<f64>::from_vecs( 5, 5, val, row_index, col_start )
}

#[test]
fn create_from_vecs_sparse2() {
    let s = test_sparse_matrix();
    assert_eq!( s.rows, 5 );
    assert_eq!( s.cols, 5 );
    assert_eq!( s.nonzero, 9 );
    assert_eq!( s.get( 0, 0 ), Some( 3.0 ) );
    assert_eq!( s.get( 1, 1 ), Some( 4.0 ) );
    assert_eq!( s.get( 2, 1 ), Some( 7.0 ) );
    assert_eq!( s.get( 0, 2 ), Some( 1.0 ) );
    assert_eq!( s.get( 2, 2 ), Some( 5.0 ) );
    assert_eq!( s.get( 0, 3 ), Some( 2.0 ) );
    assert_eq!( s.get( 2, 3 ), Some( 9.0 ) );
    assert_eq!( s.get( 4, 3 ), Some( 6.0 ) );
    assert_eq!( s.get( 4, 4 ), Some( 5.0 ) );
    assert_eq!( s.get( 3, 3 ), None );
}

#[test]
fn col_index_sparse2() {
    let s = test_sparse_matrix();
    let col_index = s.col_index();
    assert_eq!( col_index.vec, vec![ 0, 1, 1, 2, 2, 3, 3, 3, 4 ] );
}

#[test]
fn col_start_from_index_sparse2() {
    let s = test_sparse_matrix();
    let col_index = Vector::<usize>::create( vec![ 0, 1, 1, 2, 2, 3, 3, 3, 4 ] );
    let col_start = s.col_start_from_index( &col_index );
    assert_eq!( col_start, vec![ 0, 1, 3, 5, 8, 9 ] );
}

#[test]
fn scale_sparse2() {
    let mut s = test_sparse_matrix();
    s.scale( &2.0 );
    assert_eq!( s.get( 0, 0 ), Some( 6.0 ) );
    assert_eq!( s.get( 1, 1 ), Some( 8.0 ) );
    assert_eq!( s.get( 2, 1 ), Some( 14.0 ) );
    assert_eq!( s.get( 0, 2 ), Some( 2.0 ) );
    assert_eq!( s.get( 2, 2 ), Some( 10.0 ) );
    assert_eq!( s.get( 0, 3 ), Some( 4.0 ) );
    assert_eq!( s.get( 2, 3 ), Some( 18.0 ) );
    assert_eq!( s.get( 4, 3 ), Some( 12.0 ) );
    assert_eq!( s.get( 4, 4 ), Some( 10.0 ) );
    assert_eq!( s.get( 3, 3 ), None );
}

#[test]
fn multiply_sparse2() {
    let s = test_sparse_matrix();
    let v = Vector::<f64>::create( vec![ 1.0, 1.0, 1.0, 1.0, 1.0 ] );
    let result = s.multiply( &v );
    assert_eq!( result.vec, vec![ 6.0, 4.0, 21.0, 0.0, 11.0 ] );
}

#[test]
fn transpose_multiply_sparse2() {
    let s = test_sparse_matrix();
    let v = Vector::<f64>::create( vec![ 1.0, 1.0, 1.0, 1.0, 1.0 ] );
    let result = s.transpose_multiply( &v );
    assert_eq!( result.vec, vec![ 3.0, 11.0, 6.0, 17.0, 5.0 ] );
}

#[test]
fn transpose_sparse2() {
    let s = test_sparse_matrix();
    let t = s.transpose();
    assert_eq!( t.rows, 5 );
    assert_eq!( t.cols, 5 );
    assert_eq!( t.nonzero, 9 );
    assert_eq!( t.get( 0, 0 ), Some( 3.0 ) );
    assert_eq!( t.get( 1, 1 ), Some( 4.0 ) );
    assert_eq!( t.get( 1, 2 ), Some( 7.0 ) );
    assert_eq!( t.get( 2, 0 ), Some( 1.0 ) );
    assert_eq!( t.get( 2, 2 ), Some( 5.0 ) );
    assert_eq!( t.get( 3, 0 ), Some( 2.0 ) );
    assert_eq!( t.get( 3, 2 ), Some( 9.0 ) );
    assert_eq!( t.get( 3, 4 ), Some( 6.0 ) );
    assert_eq!( t.get( 4, 4 ), Some( 5.0 ) );
    assert_eq!( t.get( 3, 3 ), None );
    assert_eq!( t.get( 0, 2 ), None );
}

#[test]
fn triplet_sparse2() {
    let mut triplets = vec![];
    triplets.push( (0, 3, 2.0) );
    triplets.push( (2, 3, 9.0) );
    triplets.push( (4, 3, 6.0) );
    triplets.push( (4, 4, 5.0) );
    triplets.push( (0, 0, 3.0) );
    triplets.push( (1, 1, 4.0) );
    triplets.push( (2, 1, 7.0) );
    triplets.push( (0, 2, 1.0) );
    triplets.push( (2, 2, 5.0) );
    triplets.push( (3, 2, 1.0) );
    let s = Sparse2::<f64>::from_triplets( 5, 5, &mut triplets );
    assert_eq!( s.get( 0, 0 ), Some( 3.0 ) );
    assert_eq!( s.get( 1, 1 ), Some( 4.0 ) );
    assert_eq!( s.get( 2, 1 ), Some( 7.0 ) );
    assert_eq!( s.get( 0, 2 ), Some( 1.0 ) );
    assert_eq!( s.get( 2, 2 ), Some( 5.0 ) );
    assert_eq!( s.get( 4, 3 ), Some( 6.0 ) );
    assert_eq!( s.get( 2, 3 ), Some( 9.0 ) );
    assert_eq!( s.get( 0, 3 ), Some( 2.0 ) );
    assert_eq!( s.get( 4, 4 ), Some( 5.0 ) );
    assert_eq!( s.get( 3, 2 ), Some( 1.0 ) );
    assert_eq!( s.get( 3, 3 ), None );
}


fn test_sparse_matrix_2() -> Sparse2::<f64> {
    let mut triplets = vec![
        (0, 0, 3.0), (0, 2, 1.0), (0, 3, 2.0),
        (1, 1, 4.0),
        (2, 1, 7.0), (2, 2, 5.0), (2, 3, 9.0),
        (3, 2, 1.0),
        (4, 3, 6.0), (4, 4, 5.0)
    ];
    Sparse2::<f64>::from_triplets( 5, 5, &mut triplets )
}

#[test]
fn solve_bicg_sparse2() {
    let s = test_sparse_matrix_2();
    let b = Vector::<f64>::create( vec![ 1.0, 2.0, 3.0, 4.0, 5.0 ] );
    let mut x = Vector::<f64>::new( 5, 0.0 );
    let max_iter = 1000;
    let tol = 1e-8;
    let result = s.solve_bicg( &b, &mut x, max_iter, tol, 1 );
    assert!( result.is_ok() && !result.is_err() );
    assert!( (x[0] - 0.518518519).abs() < tol );
    assert!( (x[1] - 0.5).abs() < tol );
    assert!( (x[2] - 4.0).abs() < tol );
    assert!( (x[3] - -2.277777778).abs() < tol );
    assert!( (x[4] - 3.733333333).abs() < tol );
}

#[test]
fn solve_bicgstab_sparse2() {
    let s = test_sparse_matrix_2();
    let b = Vector::<f64>::create( vec![ 1.0, 2.0, 3.0, 4.0, 5.0 ] );
    let mut x = Vector::<f64>::new( 5, 0.0 );
    let max_iter = 1000;
    let tol = 1e-8;
    let result = s.solve_bicgstab( &b, &mut x, max_iter, tol );
    assert!( result.is_ok() && !result.is_err() );
    assert!( (x[0] - 0.518518519).abs() < tol );
    assert!( (x[1] - 0.5).abs() < tol );
    assert!( (x[2] - 4.0).abs() < tol );
    assert!( (x[3] - -2.277777778).abs() < tol );
    assert!( (x[4] - 3.733333333).abs() < tol );
}

#[test]
fn solve_cg_sparse2() {
    let mut triplets = vec![
        (0, 0, 4.0), (0, 1, 1.0), 
        (1, 0, 1.0), (1, 1, 3.0)
    ];
    let s = Sparse2::<f64>::from_triplets( 2, 2, &mut triplets );
    let b = Vector::<f64>::create( vec![ 1.0, 2.0 ] );
    let mut x = Vector::<f64>::create(vec![ 2.0, 1.0 ]);
    let max_iter = 1000;
    let tol = 1e-8;
    let result = s.solve_cg( &b, &mut x, max_iter, tol );
    assert!( result.is_ok() && !result.is_err() );
    assert!( (x[0] - ( 1. / 11. )).abs() < tol );
    assert!( (x[1] - ( 7. / 11. )).abs() < tol );
}

#[test]
fn solve_qmr_sparse2() {
    let mut triplets = vec![
        (0, 0, 4.0), (0, 1, 1.0), 
        (1, 0, 1.0), (1, 1, 3.0)
    ];
    let s = Sparse2::<f64>::from_triplets( 2, 2, &mut triplets );
    let b = Vector::<f64>::create( vec![ 1.0, 2.0 ] );
    let mut x = Vector::<f64>::create(vec![ 2.0, 1.0 ]);
    let max_iter = 1000;
    let tol = 1e-8;
    let result = s.solve_qmr( &b, &mut x, max_iter, tol );
    assert!( result.is_ok() && !result.is_err() );
    assert!( (x[0] - ( 1. / 11. )).abs() < tol );
    assert!( (x[1] - ( 7. / 11. )).abs() < tol );
}

#[test]
fn convert_to_triplets() {
    let s = test_sparse_matrix_2();
    let mut triplets = s.to_triplets();
    assert_eq!( triplets.len(), 10 );
    triplets.sort_by_key( |triplet| triplet.1 );
    let mut original = vec![
        (0, 0, 3.0), (0, 2, 1.0), (0, 3, 2.0),
        (1, 1, 4.0),
        (2, 1, 7.0), (2, 2, 5.0), (2, 3, 9.0),
        (3, 2, 1.0),
        (4, 3, 6.0), (4, 4, 5.0)
    ];
    original.sort_by_key( |triplet| triplet.1 );
    assert_eq!( triplets, original );
}

#[test]
fn insert_sparse2() {
    let mut s = test_sparse_matrix_2();
    s.insert( 3, 3, 4.0 );
    assert_eq!( s.get( 3, 3 ), Some( 4.0 ) );
    s.insert( 3, 3, 1.0 );
    assert_eq!( s.get( 3, 3 ), Some( 1.0 ) );
}

#[test]
fn convert_to_dense() {
    let sparse = test_sparse_matrix_2();
    let dense = sparse.to_dense();
    assert_eq!( dense[(0,0)], 3.0 );
    assert_eq!( dense[(0,1)], 0.0 );
    assert_eq!( dense[(0,2)], 1.0 );
    assert_eq!( dense[(0,3)], 2.0 );
    assert_eq!( dense[(0,4)], 0.0 );
    assert_eq!( dense[(1,0)], 0.0 );
    assert_eq!( dense[(1,1)], 4.0 );
    assert_eq!( dense[(1,2)], 0.0 );
    assert_eq!( dense[(1,3)], 0.0 );
    assert_eq!( dense[(1,4)], 0.0 );
    assert_eq!( dense[(2,0)], 0.0 );
    assert_eq!( dense[(2,1)], 7.0 );
    assert_eq!( dense[(2,2)], 5.0 );
    assert_eq!( dense[(2,3)], 9.0 );
    assert_eq!( dense[(2,4)], 0.0 );
    assert_eq!( dense[(3,0)], 0.0 );
    assert_eq!( dense[(3,1)], 0.0 );
    assert_eq!( dense[(3,2)], 1.0 );
    assert_eq!( dense[(3,3)], 0.0 );
    assert_eq!( dense[(3,4)], 0.0 );
    assert_eq!( dense[(4,0)], 0.0 );
    assert_eq!( dense[(4,1)], 0.0 );
    assert_eq!( dense[(4,2)], 0.0 );
    assert_eq!( dense[(4,3)], 6.0 );
    assert_eq!( dense[(4,4)], 5.0 );
}
