use ohsl::matrix::{Sparse, Sparse64, Triplet, Tr64};
use ohsl::vector::{Vec64};

#[test]
    fn test_sparse_unspecified_size() {
        let s = Sparse::<i32>::empty( 10, 10 );
        assert_eq!( s.rows(), 10 );
        assert_eq!( s.cols(), 10 );
        assert_eq!( s.nonzero(), 0 );
    }

    #[test]
    fn test_triplet() {
        let trip = Triplet::<f64>::empty();
        assert_eq!( trip.row(), 0 );
        assert_eq!( trip.col(), 0 );
        assert_eq!( trip.val(), 0.0 );
        let trip2 = Triplet::<f64>::new( 1, 2, 1.0 );
        assert_eq!( trip2.row(), 1 );
        assert_eq!( trip2.col(), 2 );
        assert_eq!( trip2.val(), 1.0 );
        assert!( trip.compare( &trip2 ));
        assert!( !trip2.compare( &trip ));
    }

    #[test]
    fn test_sparse_triplet() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        triplets.push( Tr64::new( 2, 2, 7.0 ) );
        let sparse = Sparse::<f64>::new( 3, 3, &mut triplets );
        assert_eq!( sparse.rows(), 3 ); 
        assert_eq!( sparse.cols(), 3 ); 
        assert_eq!( sparse.nonzero(), 7 ); 
        let values = sparse.val();
        assert_eq!( values[0], 1.0 );
        assert_eq!( values[1], 3.0 );
        assert_eq!( values[2], 2.0 );
        assert_eq!( values[3], 4.0 );
        assert_eq!( values[4], 6.0 );
        assert_eq!( values[5], 5.0 );
        assert_eq!( values[6], 7.0 );
        let index = sparse.row_index();
        assert_eq!( index[0], 0 );
        assert_eq!( index[1], 1 );
        assert_eq!( index[2], 0 );
        assert_eq!( index[3], 1 );
        assert_eq!( index[4], 2 );
        assert_eq!( index[5], 1 );
        assert_eq!( index[6], 2 );
        let start = sparse.col_start();
        assert_eq!( start[0], 0 );
        assert_eq!( start[1], 2 );
        assert_eq!( start[2], 5 );
        assert_eq!( start[3], 7 );
    }

    #[test]
    fn test_sparse_scale() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        triplets.push( Tr64::new( 2, 2, 7.0 ) );
        let mut sparse = Sparse::<f64>::new( 3, 3, &mut triplets );
        sparse.scale( 2.0 );
        let values = sparse.val();
        assert_eq!( values[0], 2.0 );
        assert_eq!( values[1], 6.0 );
        assert_eq!( values[2], 4.0 );
        assert_eq!( values[3], 8.0 );
        assert_eq!( values[4], 12.0 );
        assert_eq!( values[5], 10.0 );
        assert_eq!( values[6], 14.0 );
    }

    #[test]
    fn test_sparse_vector_mult() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        triplets.push( Tr64::new( 2, 2, 7.0 ) );
        let sparse = Sparse::<f64>::new( 3, 3, &mut triplets );
        let v = Vec64::new( 3, 1.0 );
        let mul = sparse.multiply( &v );
        assert_eq!( mul[0], 3.0 );
        assert_eq!( mul[1], 12.0 );
        assert_eq!( mul[2], 13.0 );
        let trans_mul = sparse.transpose_multiply( &v );
        assert_eq!( trans_mul[0], 4.0 );
        assert_eq!( trans_mul[1], 12.0 );
        assert_eq!( trans_mul[2], 12.0 );
    }

    #[test]
    fn test_sparse_solve_bicg() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        triplets.push( Tr64::new( 2, 2, 7.0 ) );
        let sparse = Sparse::<f64>::new( 3, 3, &mut triplets );
        let b = Vec64::new( 3, 1.0 );
        // Set the initial guess
        let mut x = Vec64::new( 3, 1.0 );
        x[0] = 0.136;
        x[1] = 0.432;
        x[2] = -0.227;
        let max_iter = 100;
        let tol = 1.0e-6;
        let (solution, iter) = sparse.solve_bicg( b, x, max_iter, tol);
        assert!( (solution[0] - 6.0/44.0).abs() < tol );
        assert!( (solution[1] - 19.0/44.0).abs() < tol );
        assert!( (solution[2] + 10.0/44.0).abs() < tol );
        assert!( iter < max_iter );
    }

    #[test]
    fn test_sparse_solve_bicgstab() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        triplets.push( Tr64::new( 2, 2, 7.0 ) );
        let sparse = Sparse::<f64>::new( 3, 3, &mut triplets );
        let b = Vec64::new( 3, 1.0 );
        // Set the initial guess
        let mut x = Vec64::new( 3, 1.0 );
        x[0] = 0.136;
        x[1] = 0.432;
        x[2] = -0.227;
        let max_iter = 100;
        let tol = 1.0e-6;
        let solution = sparse.solve_bicgstab( b, x, max_iter, tol);
        assert!( (solution[0] - 6.0/44.0).abs() < tol );
        assert!( (solution[1] - 19.0/44.0).abs() < tol );
        assert!( (solution[2] + 10.0/44.0).abs() < tol );
    }

    #[test]
    fn test_sparse_insert() {
        let mut triplets = Vec::new();
        triplets.push( Tr64::new( 0, 0, 1.0 ) );
        triplets.push( Tr64::new( 0, 1, 2.0 ) );
        triplets.push( Tr64::new( 1, 0, 3.0 ) );
        triplets.push( Tr64::new( 1, 1, 4.0 ) );
        triplets.push( Tr64::new( 1, 2, 5.0 ) );
        triplets.push( Tr64::new( 2, 1, 6.0 ) );
        let mut sparse = Sparse64::new( 3, 3, &mut triplets );
        sparse.insert( 0, 0, 1.0 );
        sparse.insert( 2, 2, 7.0 );
        let col_index = sparse.col_index();
        assert_eq!( col_index[0], 0 );
        assert_eq!( col_index[1], 0 );
        assert_eq!( col_index[2], 1 );
        assert_eq!( sparse.val()[4], 6.0 );
        assert_eq!( sparse.val()[5], 5.0 );
        assert_eq!( sparse.val()[6], 7.0 );
        let b = Vec64::new( 3, 1.0 );
        let mut x = Vec64::new( 3, 1.0 );
        x[0] = 0.136;
        x[1] = 0.432;
        x[2] = -0.227;
        let max_iter = 100;
        let tol = 1.0e-6;
        let solution = sparse.solve_bicgstab( b, x, max_iter, tol);
        assert!( (solution[0] - 6.0/44.0).abs() < tol );
        assert!( (solution[1] - 19.0/44.0).abs() < tol );
        assert!( (solution[2] + 10.0/44.0).abs() < tol );
    }
