use std::mem;
use core::ops::{Index, IndexMut};
pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};
pub use crate::traits::{Number, Signed, Zero, One};

impl<T> Index<usize> for Matrix<T> {
    type Output = Vector<T>;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: usize ) -> &'a Vector<T> {
        &self.mat[ index ]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: usize ) -> &mut Vector<T> {
        &mut self.mat[ index ] 
    }
}

impl<T> Matrix<T> {
    /// Remove all the elements from the matrix 
    #[inline]
    pub fn clear(&mut self) {
        self.mat.clear();
        self.rows = 0;
        self.cols = 0;
    }
}

impl<T: Clone + Copy + Number> Matrix<T> {
    /// Get a row of the matrix as a vector
    #[inline]
    pub fn get_row(&self, row: usize ) -> Vector<T> {
        if self.rows <= row { panic!( "Matrix range error in get_row" ); }
        self.mat[ row ].clone()
    }

    /// Get a column of the matrix as a vector 
    #[inline]
    pub fn get_col(&self, col: usize ) -> Vector<T> {
        if self.rows <= col { panic!( "Matrix range error in get_col" ); }
        let mut result = Vector::<T>::new( self.rows, T::zero() );
        for i in 0..self.rows {
            result[ i ] = self[i][col]
        }
        result
    }

    /// Set a row of the matrix using a vector 
    #[inline]
    pub fn set_row(&mut self, row: usize, vec: Vector<T> ) {
        if vec.size() != self.cols { panic!( "Matrix size error in set_row" ); }
        if self.rows <= row { panic!( "Matrix range error in set_row" ); }
        self[ row ] = vec;
    }

    /// Set a column of the matrix using a vector 
    #[inline]
    pub fn set_col(&mut self, col: usize, vec: Vector<T> ) {
        if vec.size() != self.rows { panic!( "Matrix size error in set_col" ); }
        if self.rows <= col { panic!( "Matrix range error in set_col" ); }
        for i in 0..self.rows {
            self[i][col] = vec[i];
        }
    }

    /// Delete a row from the matrix 
    #[inline]
    pub fn delete_row(&mut self, row: usize ) {
        if self.rows <= row { panic!( "Matrix range error in delete_row" ); }
        self.mat.remove( row );
        self.rows -= 1;
    }

    /// Multiply the matrix by a (column) vector and return a vector 
    #[inline]
    pub fn multiply(&self, vec: &Vector<T> ) -> Vector<T> {
        if vec.size() != self.cols { panic!( "Matrix dimensions do not agree in multiply." ); }
        let mut result = Vector::<T>::empty();
        for row in 0..self.rows {
           result.push( self[row].dot( vec ) );
        }
        result
    }

    /// Create a square identity matrix of specified size 
    #[inline]
    pub fn eye( size: usize ) -> Self {
        let mut identity = Matrix::<T>::new( size, size, T::zero() ); 
        for i in 0..size {
            identity[i][i] = T::one();
        }
        identity
    }

    /// Resize the matrix (empty entries are appended if necessary)
    #[inline]
    pub fn resize(&mut self, n_rows: usize, n_cols: usize ) {
        let temp = self.clone();
        *self = Matrix::<T>::new( n_rows, n_cols, T::zero() );
        for i in 0..n_rows {
            for j in 0..n_cols {
                if i < temp.rows() && j < temp.cols() {
                    self[i][j] = temp[i][j].clone();
                }
            }
        }
    }

    /// Transpose the matrix in place 
    #[inline]
    pub fn transpose_in_place(&mut self) {
        if self.rows == self.cols {
            for i in 0..self.rows {
                for j in i+1..self.cols {
                    let mut temp = self[i][j].clone();
                    mem::swap( &mut self[j][i], &mut temp );
                    self[i][j] = temp;
                }
            }
        } else {
            let row = Vector::<T>::new( self.rows, T::zero() );
            let mut temp = Vec::new();
            for _i in 0..self.cols {
                temp.push( row.clone() );
            }
            for i in 0..self.rows {
                for j in 0..self.cols {
                    temp[j][i] = self[i][j].clone();
                }
            }
            self.mat = temp;
            mem::swap( &mut self.rows, &mut self.cols );
        }
    }

    /// Return the transpose of the matrix 
    #[inline]
    pub fn transpose(&self) -> Matrix<T> {
        let mut temp: Matrix<T> = self.clone();
        temp.transpose_in_place();
        temp
    }

    /// Swap two rows of the matrix 
    #[inline]
    pub fn swap_rows(&mut self, row_1: usize, row_2: usize ) {
        if self.rows <= row_1 || self.rows <= row_2 { panic!( "Matrix swap row range error." ); } 
        let mut temp = self.mat[ row_1 ].clone();
        mem::swap( &mut self.mat[ row_2 ], &mut temp );
        self.mat[ row_1 ] = temp;
    }

    /// Swap two elements of the matrix 
    #[inline]
    pub fn swap_elem(&mut self, row_1: usize, col_1: usize, row_2: usize, col_2: usize ) {
        let mut temp = self[row_1][col_1].clone();
        mem::swap( &mut self[row_2][col_2], &mut temp );
        self[row_1][col_1] = temp;
    }

    /// Fill the matrix with specified elements
    #[inline]
    pub fn fill(&mut self, elem: T ) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] = elem.clone();
            }
        }
    }

    /// Fill the leading diagonal of the matrix with specified elements
    #[inline]
    pub fn fill_diag(&mut self, elem: T ) {
        let n: usize = if self.cols < self.rows { self.cols } else { self.rows };
        for i in 0..n {
            self[i][i] = elem.clone();
        }
        
    }

    /// Fill a diagonal band of the matrix with specified elements
    /// ( offset above main diagonal +, below main diagonal - )
    #[inline]
    pub fn fill_band(&mut self, offset: isize, elem: T ) {
        for row in 0..self.rows {
            let i = (row as isize) + offset; 
            if (i as usize) < self.cols &&  i >= 0 {
                self[ row ][ i as usize ] = elem.clone();
            } 
        }
    }

    /// Fill the main three diagonals of the matrix with specified elements 
    #[inline]
    pub fn fill_tridiag(&mut self, lower: T, diag: T, upper: T ) {
        self.fill_band( -1, lower );
        self.fill_diag( diag );
        self.fill_band( 1, upper );
    }

    /// Fill a row of the matrix with specified elements
    #[inline]
    pub fn fill_row(&mut self, row: usize, elem: T ) {
        if self.rows <= row { panic!( "Matrix range error in fill_row" ); }
        for j in 0..self.cols {
            self[row][j] = elem.clone();
        }
    }

    /// Fill a column of the matrix with specified elements
    #[inline]
    pub fn fill_col(&mut self, col: usize, elem: T ) {
        if self.cols <= col { panic!( "Matrix range error in fill_col" ); }
        for i in 0..self.rows {
            self[i][col] = elem.clone();
        }
    }
}