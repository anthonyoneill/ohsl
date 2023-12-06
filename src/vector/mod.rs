pub mod arithmetic;
pub mod functions;
pub mod operations;
pub mod vec_f64;
pub mod vec_cmplx;

use std::{fmt, fs::File, io::Write};
pub use crate::traits::{Number, Zero, One};
pub use crate::complex::Complex;

#[derive(PartialEq)]
pub struct Vector<T> {
    pub vec: Vec<T>,
    size: usize,
}

pub type Vec64 = Vector<f64>;

impl<T> Vector<T> {
    /// Create a new vector of unspecified size
    #[inline]
    pub const fn empty() -> Self {
        let vec = Vec::<T>::new();
        let size = 0;
        Vector { vec, size }
    }

    /// Create a vector from an std::vec::Vec
    #[inline]
    pub fn create( vec: Vec<T> ) -> Self {
        let size = vec.len();
        Vector { vec, size }
    }

    /// Return the size of the vector 
    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }
}

impl<T: Clone> Vector<T> {
    /// Create a new vector of specified size
    #[inline]
    pub fn new( size: usize, elem: T ) -> Self {
        let vec = vec![ elem; size ];
        Vector { vec, size }
    }
}

impl<T: Clone + Number> Vector<T> {
    /// Create a vector of zeros 
    #[inline]
    pub fn zeros( size: usize ) -> Self {
        let vec = vec![ T::zero(); size ];
        Vector{ vec, size }
    }

    /// Create a vector of ones
    #[inline]
    pub fn ones( size: usize ) -> Self {
        let vec = vec![ T::one(); size ];
        Vector{ vec, size }
    }
}

impl<T: Clone> Clone for Vector<T> {
    /// Clone the vector
    #[inline]
    fn clone(&self) -> Self {
        Self::create( self.vec.clone() )
    }
}

impl<T: fmt::Display> Vector<T> {
    /// Print the vector to a file
    #[inline]
    pub fn output(&self, filename: &str) {
        let mut f = File::create(filename).expect("Unable to create file");
        for i in 0..self.size {                                                                                                                                                                  
            writeln!(f, "{}", self.vec[i]).unwrap();                                                                                                                            
        }
    }
}

impl<T> fmt::Debug for Vector<T> where
    T: fmt::Debug
{
    /// Format the output [ v_0, v_1, ... ]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vec )
    }
}

impl<T> fmt::Display for Vector<T> where
    T: fmt::Debug
{
    /// Format the output [ v_0, v_1, ... ]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vec )
    }
} 