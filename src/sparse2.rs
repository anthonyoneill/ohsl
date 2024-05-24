use crate::vector::Vector;
pub use crate::traits::{Number, Zero};

pub struct Sparse2<T> {
    pub rows: usize,                // Number of rows
    pub cols: usize,                // Number of columns
    pub nonzero: usize,             // Number of non-zero elements
    pub val: Vec<T>,                // Non-zero values
    pub row_index: Vec<usize>,      // Row indices of non-zero values
    pub col_start: Vec<usize>,      // Column indices of non-zero values
}

impl<T: Copy + Number> Sparse2<T> {
    /// Create a new sparse matrix of specified size
    pub fn new( rows: usize, cols: usize ) -> Self {
        Self {
            rows,
            cols,
            nonzero: 0,
            //val: Vec::new(),
            val: vec![ T::zero(); 0 ],
            //row_index: Vec::new(),
            row_index: vec![ 0; 0 ],
            col_start: vec![ 0; cols + 1],
        }
    }

    // Create a new sparse matrix of specified size and number of non-zero elements
    fn new_nonzero( rows: usize, cols: usize, nonzero: usize ) -> Self {
        Self {
            rows,
            cols,
            nonzero,
            val: vec![ T::zero(); nonzero ],
            row_index: vec![ 0; nonzero ],
            col_start: vec![ 0; cols + 1 ],
        }
    } 

    /// Create a new sparse matrix of specified size using a vector of values
    pub fn from_vecs( rows: usize, cols: usize, val: Vec<T>, row_index: Vec<usize>, col_start: Vec<usize> ) -> Self {
        //TODO check that the vectors are the correct length val.len() == row_index.len()
        Self {
            rows: rows,
            cols: cols,
            //nonzero: val.len(),
            nonzero: col_start[ col_start.len() - 1 ],
            val,
            row_index,
            col_start,
        }
    }

    /// Return a vector of column indices relating to each value ( triplet form )
    pub fn col_index( &self ) -> Vector<usize> {
        let mut temp = Vector::empty();
        if self.nonzero == 0 { return temp; }
        if self.col_start.len() < self.cols + 1 { 
            panic!( "Sparse matrix col_index: Some columns have no entries." ); 
        }
        let mut gaps = vec![ 0; self.col_start.len() - 1 ];
        for k in 0..gaps.len() {
            gaps[ k ] = self.col_start[ k + 1 ] - self.col_start[ k ];
            for _j in 0..gaps[ k ] {
                temp.push( k );
            }
        }
        temp
    }

    /// Get a specific element from the sparse matrix if it exists
    pub fn get( &self, row: usize, col: usize ) -> Option<T> {
        if self.rows <= row { panic!( "Sparse matrix get: row range error." ); }
        if self.cols <= col { panic!( "Sparse matrix get: col range error." ); }
        if self.col_start.len() <= col { 
            panic!( "Sparse matrix get: Some columns have no entries." ); 
        }
        // Convert to triplet format
        let col_index = self.col_index();
        // Check if the element exists
        for k in 0..self.nonzero {
            if ( self.row_index[ k ] == row ) && ( col_index[ k ] == col ) {
                return Some( self.val[ k ] );
            }
        }
        None
    }
    //TODO maybe index operator should be implemented 

    /// Calculate column start vector for a given column index vector
    pub fn col_start_from_index( &self, col_index: &Vector<usize> ) -> Vec<usize> {
        let mut col_start = vec![ 0; self.cols + 1 ];
        // Compute the number of non-zero elements in each column
        for n in 0..self.nonzero {
            col_start[ col_index[ n ] ] += 1;
        }
        let mut sum = 0;
        // Cumulative sum to get the start of each column
        for k in 0..self.cols {
            let ck = col_start[ k ];
            col_start[ k ] = sum;
            sum += ck;
        }
        col_start[ self.cols ] = sum;
        col_start
    }

    /// Insert a non-zero element into the sparse matrix at the specified row and column
    pub fn insert( &mut self, row: &usize, col: &usize, value: &T ) {
        if self.rows <= *row { panic!( "Sparse matrix insert: row range error." ); }
        if self.cols <= *col { panic!( "Sparse matrix insert: col range error." ); }
        if self.col_start.len() <= *col { 
            panic!( "Sparse matrix insert: Some columns have no entries." ); 
        }
        let mut col_index = self.col_index();
        // Check if the element already exists
        for k in 0..self.nonzero {
            if ( self.row_index[ k ] == *row ) && ( col_index[ k ] == *col ) {
                self.val[ k ] = *value;
                return;
            }
        }
        if self.nonzero == 0 {
            self.row_index.push( *row );
            self.val.push( *value );
            col_index.push( *col );
            self.nonzero += 1;
            self.col_start = self.col_start_from_index( &col_index );
            return;
        }
        // If nonzero != 0 insert new element
        let find_col = col_index.find( *col );
        col_index.insert( find_col, *col );
        let mut find_row = find_col;
        while self.row_index[ find_row ] < *row {
            find_row += 1;
        }
        self.row_index.insert( find_row, *row );
        self.val.insert( find_row, *value );
        // Convert back to compressed sparse column format
        self.col_start = self.col_start_from_index( &col_index );
        self.nonzero += 1;
    }

    /// Scale the non-zero elements of sparse matrix by a given value
    pub fn scale( &mut self, value: &T ) {
        for k in 0..self.nonzero {
            self.val[ k ] *= *value;
        }
    }

    /// Multiply the sparse matrix A by a Vector x to the right and return the result vector A * x
    pub fn multiply( &self, x: &Vector<T> ) -> Vector<T> {
        if self.cols != x.size() { 
            panic!( "Sparse matrix multiply: matrix and vector sizes do not agree." ); 
        }
        let mut result = Vector::create( vec![ T::zero(); self.rows ] );
        for j in 0..self.cols {
            let xj = x[ j ];
            for k in self.col_start[ j ]..self.col_start[ j + 1 ] {
                result[ self.row_index[ k ] ] += self.val[ k ] * xj;
            }
        }
        result
    }

    /// Multiply the transpose matrix A^T by a Vector x and return the result vector A^T * x
    pub fn transpose_multiply( &self, x: &Vector<T> ) -> Vector<T> {
        if self.rows != x.size() { 
            panic!( "Sparse matrix transpose_multiply: matrix and vector sizes do not agree." ); 
        }
        let mut result = Vector::create( vec![ T::zero(); self.cols ] );
        for i in 0..self.cols {
            for k in self.col_start[ i ]..self.col_start[ i + 1 ] {
                result[ i ] += self.val[ k ] * x[ self.row_index[ k ] ];
            }
            
        }
        result
    }

    /// Return the transpose of the sparse matrix
    pub fn transpose( &self ) -> Sparse2<T> {
        let mut at = Sparse2::new_nonzero( self.cols, self.rows, self.nonzero );
        let mut count = vec![ 0; self.rows ];
        for i in 0..self.cols {
            for j in self.col_start[ i ]..self.col_start[ i + 1 ] {
                count[ self.row_index[ j ] ] += 1;
            }
        }
        for j in 0..self.rows {
            at.col_start[ j + 1 ] = at.col_start[ j ] + count[ j ];
        }
        count = vec![ 0; self.rows ];
        for i in 0..self.cols {
            for j in self.col_start[ i ]..self.col_start[ i + 1 ] {
                let k = self.row_index[ j ];
                let index = at.col_start[ k ] + count[ k ];
                at.row_index[ index ] = i;
                at.val[ index ] = self.val[ j ];
                count[ k ] += 1;
            }
        }
        at
    }

    //TODO solve_bicg, diagonal_preconditioner, identity_preconditioner 
    // solve_bicgstab, solve_cg, solve_qmr, solve_gmres
}