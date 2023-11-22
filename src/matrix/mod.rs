pub mod operations;
pub mod arithmetic;
pub mod solve;
pub mod functions;

use core::ops::{Index, IndexMut};
use std::{fmt, fs::File, io::Write};

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Complex;
pub use crate::vector::{Vector, Vec64};

pub struct Matrix<T> {
    mat: Vec< Vector<T> >,
    rows: usize,
    cols: usize,
}

pub type Mat64 = Matrix<f64>;

impl<T> Matrix<T> {
    /// Create a new matrix of unspecified size
    #[inline]
    pub fn empty() -> Self {
        let row = Vector::<T>::empty();
        let mut mat = Vec::new();
        mat.push( row );
        Matrix { mat, rows: 0, cols: 0 }
    }

    /// Create a matrix from an `std::vec::Vec<Vector<T>>`
    #[inline]
    pub fn create( mat: Vec< Vector<T> > ) -> Self {
        let rows = mat.len();
        let cols = mat[0].size();
        Matrix { mat, rows, cols }
    }

    /// Return the number of rows in the matrix 
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Return the number of columns in the matrix 
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Return the number of elements in the matrix 
    #[inline]
    pub fn numel(&self) -> usize {
        self.cols * self.rows
    }
}

impl<T: Clone + Number> Matrix<T> {
    /// Create a new matrix of specified size
    #[inline]
    pub fn new( rows: usize, cols: usize, elem: T ) -> Self {
        let row = Vector::<T>::new( cols, elem );
        let mut mat = Vec::new();
        for _i in 0..rows {
            mat.push( row.clone() );
        }
        Matrix { mat, rows, cols }
    }
}

impl<T: Clone> Clone for Matrix<T> {
    /// Clone the matrix
    #[inline]
    fn clone(&self) -> Self {
        Self::create( self.mat.clone() )
    }
}

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

impl<T> fmt::Debug for Matrix<T> where
    T: fmt::Debug
{
    /// Format the output 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows-1 {
            writeln!(f, "\t{:?}", self.mat[i] ).unwrap();
        }
        write!(f, "\t{:?}", self.mat[self.rows-1] )
    }
}

impl<T> fmt::Display for Matrix<T> where
    T: fmt::Debug
{
    /// Format the output 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows-1 {
            writeln!(f, "\t{:?}", self.mat[i] ).unwrap();
        }
        write!(f, "\t{:?}", self.mat[self.rows-1] )
    }
} 

impl<T: fmt::Display> Matrix<T> {
    /// Print the matrix to a file
    #[inline]
    pub fn output(&self, filename: &str) {
        let mut f = File::create(filename).expect("Unable to create file");
        for i in 0..self.rows {  
            for j in 0..self.cols {
                write!(f, "\t{}", self[i][j] ).unwrap();
            }                                                                                                                                                                
            writeln!(f, "").unwrap();                                                                                                                            
        }
    }
}