use core::ops::{Index, IndexMut};

pub struct Polynomial<T> {
    /// The coefficients of the polynomial
    coeffs: Vec<T>,
}

impl<T> Polynomial<T> {
    /// Create a new polynomial of unspecified size
    #[inline]
    pub const fn empty() -> Self {
        let coeffs = Vec::<T>::new();
        Polynomial { coeffs }
    }

    /// Create a new polynomial
    #[inline]
    pub fn new( coeffs: Vec<T> ) -> Self {
        Polynomial { coeffs }
    }
    
    /// Return the size of the polynomial
    #[inline]
    pub fn size(&self) -> usize {
        self.coeffs.len()
    }
}

impl<T> Index<usize> for Polynomial<T> {
    type Output = T;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: usize ) -> &'a T {
        &self.coeffs[ index ]
    }
}

impl<T> IndexMut<usize> for Polynomial<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: usize ) -> &mut T {
        &mut self.coeffs[ index ] 
    }
}