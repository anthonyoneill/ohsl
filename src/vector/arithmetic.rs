use core::ops::{Neg, Add, Sub, Mul, Div};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use crate::traits::Number;
use crate::vector::Vector;

impl<T: Clone + Neg<Output = T>> Neg for Vector<T> {
    type Output = Self;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        let mut result = self.vec.clone();
        for i in 0..result.len() {
            result[i] = -result[i].clone();
        }
        Self::Output::create( result )
    }
}

// Non-consuming addition
impl<T: Clone + Number + Copy> Add<&Vector<T>> for Vector<T> {
    type Output = Self;
    /// Add the elements of two vectors together ( binary + )
    #[inline]
    fn add(self, plus: &Self) -> Self::Output {
        if self.size() != plus.size() { panic!( "Vector sizes do not agree (+)." ); }
        let mut result = Vec::new();
        for i in 0..self.size() {
            result.push( self.vec[i] + plus.vec[i] );
        }
        Self::Output::create( result )
    }
}

// Consuming addition
impl<T: Clone + Number + Copy> Add<Vector<T>> for Vector<T> {
    type Output = Self;
    /// Add the elements of two vectors together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        self + &plus
    }
}

// Non-consuming subtraction
impl<T: Clone + Number + Copy> Sub<&Vector<T>> for Vector<T> {
    type Output = Self;
    /// Subtract the elements of one vector from another ( binary - )
    #[inline]
    fn sub(self, minus: &Self) -> Self::Output {
        if self.size() != minus.size() { panic!( "Vector sizes do not agree (-)." ); }
        let mut result = Vec::new();
        for i in 0..self.size() {
            result.push( self.vec[i] - minus.vec[i] );
        }
        Self::Output::create( result )
    }
}

// Consuming subtraction
impl<T: Clone + Number + Copy> Sub<Vector<T>> for Vector<T> {
    type Output = Self;
    /// Subtract the elements of one vector from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        self - &minus
    }
}

impl<T: Clone + Number> Mul<T> for Vector<T> {
    type Output = Self;
    /// Multiply a vector by a scalar (vector * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..self.size() {
            result.push( self.vec[i].clone() * scalar.clone() );
        }
        Self::Output::create( result )
    }
}

impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;
    /// Allow multiplication on the left by f64 (f64 * vector)
    #[inline]
    fn mul(self, vector: Vector<f64>) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..vector.size() {
            result.push( self.clone() * vector.vec[i].clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> Div<T> for Vector<T> {
    type Output = Self;
    /// Divide a vector by a scalar (vector / scalar)
    fn div(self, scalar: T) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..self.size() {
            result.push( self.vec[i].clone() / scalar.clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> AddAssign for Vector<T> {
    /// Add a vector to a mutable vector and assign the result ( += )
    fn add_assign(&mut self, rhs: Self) {
        if self.size() != rhs.size() { panic!( "Vector sizes do not agree (+=)." ); }
        for i in 0..self.size() {
            self.vec[i] += rhs.vec[i].clone();
        }
    }
}

impl<T: Clone + Number> SubAssign for Vector<T> {
    /// Substract a vector from a mutable vector and assign the result ( -= )
    fn sub_assign(&mut self, rhs: Self) {
        if self.size() != rhs.size() { panic!( "Vector sizes do not agree (-=)." ); }
        for i in 0..self.size() {
            self.vec[i] -= rhs.vec[i].clone();
        }
    }
}

impl<T: Clone + Number> AddAssign<T> for Vector<T> {
    /// Add the same value to every element in a mutable vector 
    fn add_assign(&mut self, rhs: T) {
        for i in 0..self.size() {
            self.vec[i] += rhs.clone();
        }
    }
} 

impl<T: Clone + Number> SubAssign<T> for Vector<T> {
    /// Subtract the same value from every element in a mutable vector 
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..self.size() {
            self.vec[i] -= rhs.clone();
        }
    }
} 

impl<T: Clone + Number> MulAssign<T> for Vector<T> {
    /// Multiply every element in a mutable vector by a scalar 
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.size() {
            self.vec[i] *= rhs.clone();
        }
    }
} 

impl<T: Clone + Number> DivAssign<T> for Vector<T> {
    /// Divide every element in a mutable vector by a scalar 
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.size() {
            self.vec[i] /= rhs.clone();
        }
    }
}
