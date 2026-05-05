use core::ops::{Index, IndexMut};
pub use crate::triad::Triad;
pub use crate::matrix::{Matrix, Mat64};
pub use crate::traits::{Number, Signed, Zero, One};

impl<T> Index<(usize, usize, usize)> for Triad<T> {
    type Output = T;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: (usize, usize, usize) ) -> &'a T {
        //&self.tri[ index.0 * self.rows * self.cols + index.1 * self.cols + index.2 ]
        &self.tri[ self.cols * ( index.0 * self.rows + index.1 ) + index.2 ]
    }
}

impl<T> IndexMut<(usize, usize, usize)> for Triad<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: (usize, usize, usize) ) -> &mut T {
        //&mut self.tri[ index.0 * self.rows * self.cols + index.1 * self.cols + index.2 ]
        &mut self.tri[ self.cols * ( index.0 * self.rows + index.1 ) + index.2 ]
    }
}

impl<T> Triad<T> {
    /// Remove all the elements from the matrix 
    #[inline]
    pub fn clear(&mut self) {
        self.tri.clear();
        self.panels = 0;
        self.rows = 0;
        self.cols = 0;
    }
}

impl<T: Clone + Copy + Number> Triad<T> {
    /// Get a panel of the triad as a Matrix
    #[inline]
    pub fn get_panel(&self, panel: usize ) -> Matrix<T> {
        if self.panels <= panel { panic!( "Triad range error in get_panel" ); }
        let mut result = Matrix::<T>::new( self.rows, self.cols, T::zero() );
        let size = self.rows * self.cols;
        for j in 0..size {
            result.mat[ j ] = self.tri[ self.cols * self.rows * panel + j ];
        }
        result
    }

    /*// Multiply the matrix by a (column) vector and return a vector 
    #[inline]
    pub fn multiply(&self, vec: &Vector<T> ) -> Vector<T> {
        if vec.size() != self.cols { panic!( "Matrix dimensions do not agree in multiply." ); }
        let mut result = Vector::<T>::empty();
        for row in 0..self.rows {
           result.push( self.get_row( row ).dot( vec ) );
        }
        result
    }*/

    //TODO 
    /*// Multiply the triad by a (column) vector and return a matrix */

    //TODO Triad * vector * vector -> vector 
}