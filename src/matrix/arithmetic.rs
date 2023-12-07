use core::ops::{Neg, Add, Sub, Mul, Div};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};
pub use crate::traits::{Number, Signed, Zero, One};

// Non-consuming negation
impl<T: Copy + Neg<Output = T> + Signed> Neg for &Matrix<T> {
    type Output = Matrix<T>;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = -self[i][j];
            }
        }
        result
    }
}

// Consuming negation
impl<T: Copy + Neg<Output = T> + Signed> Neg for Matrix<T> {
    type Output = Self;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

// Non-consuming addition
impl<T: Copy + Number> Add<&Matrix<T>> for &Matrix<T> {
    type Output = Matrix<T>;
    /// Add the elements of two matrices together ( binary + )
    #[inline]
    fn add(self, plus: &Matrix<T>) -> Self::Output {
        if self.rows != plus.rows { panic!( "Matrix row dimensions do not agree (+)." ); }
        if self.cols != plus.cols { panic!( "Matrix col dimensions do not agree (+)." ); }
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j] + plus[i][j];
            }
        }
        result
    }
}

// Consuming addition
impl<T: Copy + Number> Add<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Add the elements of two matrices together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        &self + &plus
    }
}

// Non-consuming subtraction
impl<T: Copy + Number> Sub<&Matrix<T>> for &Matrix<T> {
    type Output = Matrix<T>;
    /// Subtract the elements of one matrix from another ( binary - )
    #[inline]
    fn sub(self, minus: &Matrix<T>) -> Self::Output {
        if self.rows != minus.rows { panic!( "Matrix row dimensions do not agree (-)." ); }
        if self.cols != minus.cols { panic!( "Matrix col dimensions do not agree (-)." ); }
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j] - minus[i][j];
            }
        }
        result
    }
}

// Consuming subtraction
impl<T: Copy + Number> Sub<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Subtract the elements of one matrix from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        &self - &minus
    }
}

// Non-consuming scalar multiplication
impl<T: Copy + Number> Mul<T> for &Matrix<T> {
    type Output = Matrix<T>;
    /// Multiply a matrix by a scalar (matrix * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j] * scalar;
            }
        }
        result
    }
}

// Consuming scalar multiplication
impl<T: Copy + Number> Mul<T> for Matrix<T> {
    type Output = Self;
    /// Multiply a matrix by a scalar (matrix * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        &self * scalar
    }
}

impl Mul<Matrix<f64>> for f64 {
    type Output = Matrix<f64>;
    /// Allow multiplication on the left by f64 (f64 * matrix)
    #[inline]
    fn mul(self, matrix: Matrix<f64>) -> Self::Output {
        let mut result = Matrix::<f64>::new( matrix.rows(), matrix.cols(), 0.0 );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = matrix[i][j] * self.clone();
            }
        }
        result
    }
}

// Non-consuming scalar division
impl<T: Copy + Number> Div<T> for &Matrix<T> {
    type Output = Matrix<T>;
    /// Divide a matrix by a scalar (matrix / scalar)
    fn div(self, scalar: T) -> Self::Output {
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j] / scalar;
            }
        }
        result
    }
}

// Consuming scalar division
impl<T: Copy + Number> Div<T> for Matrix<T> {
    type Output = Self;
    /// Divide a matrix by a scalar (matrix / scalar)
    fn div(self, scalar: T) -> Self::Output {
        &self / scalar
    }
}

// Non-consuming addition assignment
impl<T: Copy + Number> AddAssign<&Matrix<T>> for Matrix<T> {
    /// Add a matrix to a mutable matrix and assign the result ( += )
    fn add_assign(&mut self, rhs: &Self) {
        if self.rows != rhs.rows { panic!( "Matrix row dimensions do not agree (+=)." ); }
        if self.cols != rhs.cols { panic!( "Matrix col dimensions do not agree (+=)." ); }
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] += rhs[i][j];
            }
        }
    }
}

// Consuming addition assignment
impl<T: Copy + Number> AddAssign for Matrix<T> {
    /// Add a matrix to a mutable matrix and assign the result ( += )
    fn add_assign(&mut self, rhs: Self) {
        /*if self.rows != rhs.rows { panic!( "Matrix row dimensions do not agree (+=)." ); }
        if self.cols != rhs.cols { panic!( "Matrix col dimensions do not agree (+=)." ); }
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] += rhs[i][j];
            }
        }*/
        *self += &rhs
    }
}

// Non-consuming subtraction assignment
impl<T: Copy + Number> SubAssign<&Matrix<T>> for Matrix<T> {
    /// Subtract a matrix from a mutable matrix and assign the result ( -= )
    fn sub_assign(&mut self, rhs: &Self) {
        if self.rows != rhs.rows { panic!( "Matrix row dimensions do not agree (-=)." ); }
        if self.cols != rhs.cols { panic!( "Matrix col dimensions do not agree (-=)." ); }
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] -= rhs[i][j];
            }
        }
    }
}

// Consuming subtraction assignment
impl<T: Copy + Number> SubAssign for Matrix<T> {
    /// Subtract a matrix from a mutable matrix and assign the result ( -= )
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs
    }
}

impl<T: Copy + Number> MulAssign<T> for Matrix<T> {
    /// Multiply a mutable matrix by a scalar (matrix *= scalar)
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] *= rhs;
            }
        }
    }
} 

impl<T: Copy + Number> DivAssign<T> for Matrix<T> {
    /// Divide a mutable matrix by a scalar (matrix /= scalar)
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] /= rhs;
            }
        }
    }
}

impl<T: Copy + Number> AddAssign<T> for Matrix<T> {
    /// Add the same value to every element in a mutable matrix
    fn add_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] += rhs;
            }
        }
    }
}

impl<T: Copy + Number> SubAssign<T> for Matrix<T> {
    /// Subtract the same value from every element in a mutable matrix
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] -= rhs;
            }
        }
    }
}

// Consuming version of matrix-matrix multiplication
impl<T: Clone + Copy + Number> Mul<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Multiply two matrices together ( matrix * matrix )
    #[inline]
    fn mul(self, mul: Self) -> Self::Output {
        &self * &mul
    }
}

// Non-consuming version of matrix-matrix multiplication
impl<T: Clone + Copy + Number> Mul<&Matrix<T>> for &Matrix<T> {
    type Output = Matrix<T>;
    /// Multiply two matrices together ( matrix * matrix )
    #[inline]
    fn mul(self, mul: &Matrix<T> ) -> Matrix<T> {
        if self.cols != mul.rows { panic!( "Matrix dimensions do not agree (*)." ); }
        let mut result = Matrix::<T>::new( self.rows(), mul.cols(), T::zero() );
        for col in 0..mul.cols() {
            result.set_col( col, self.multiply( &mul.get_col( col ) ) );
        }
        result
    }
}

// Consuming version of matrix-vector multiplication
impl<T: Clone + Copy + Number> Mul<Vector<T>> for Matrix<T> {
    type Output = Vector<T>;
    /// Multiply a matrix with a (column) vector ( matrix * vector )
    #[inline]
    fn mul(self, vec: Vector<T> ) -> Vector<T> {
        self.multiply( &vec )
    }
}

// Non-consuming version of matrix-vector multiplication
impl<T: Clone + Copy + Number> Mul<&Vector<T>> for &Matrix<T> {
    type Output = Vector<T>;
    /// Multiply a matrix with a (column) vector ( matrix * vector )
    #[inline]
    fn mul(self, vec: &Vector<T> ) -> Vector<T> {
        self.multiply( vec )
    }
}