use ohsl::sparse_matrix::{Sparse, MatrixState, Axis};
use ohsl::vector::Vec64;

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
