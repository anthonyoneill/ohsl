pub use crate::vector::Vector;
pub use crate::traits::{Number, Signed, Zero};
pub use crate::complex::Complex;

impl<T: Clone + Signed> Vector<Complex::<T>> {
    /// Return a vector containing the complex conjugate of the elements
    #[inline]
    pub fn conj(&self) -> Self {
        let size = self.size;
        let mut vec = vec![ Complex::<T>::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().conj();
        }
        Vector { vec, size }
    }
}

impl<T: Clone + Number> Vector<Complex::<T>> {
    /// Return a vector containing the real parts of a vector of complex numbers
    #[inline]
    pub fn real(&self) -> Vector<T> {
        let size = self.size;
        let mut vec = vec![ T::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().real;
        }
        Vector { vec, size }
    } 
}