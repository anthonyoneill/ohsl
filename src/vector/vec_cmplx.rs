pub use crate::vector::Vector;
pub use crate::traits::{Number, Signed, Zero};
pub use crate::complex::Complex;

impl<T: Clone + Signed> Vector<Complex::<T>> {
    /// Return a vector containing the complex conjugate of the elements
    #[inline]
    pub fn conj(&self) -> Self {
        let size = self.size();
        let mut vec = vec![ Complex::<T>::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().conj();
        }
        Vector { vec }
    }
}

impl<T: Clone + Number> Vector<Complex::<T>> {
    /// Return a vector containing the real parts of a vector of complex numbers
    #[inline]
    pub fn real(&self) -> Vector<T> {
        let size = self.size();
        let mut vec = vec![ T::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().real;
        }
        Vector { vec }
    } 
}

impl Vector<Complex::<f64>> {
    /// Return the Inf norm: largest absolute value element (p -> infinity)
    #[inline]
    pub fn norm_inf(&self) -> f64 {
        let mut result = self.vec[0].abs();
        for i in 1..self.size() {
            if result < self.vec[i].abs() {
                result = self.vec[i].abs();
            }
        }
        result
    }
}

