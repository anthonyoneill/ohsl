use core::ops::{Add, Div, Mul, Neg, Sub};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::{fmt, cmp::Ordering};

pub use crate::traits::{Number, Signed, Zero, One, Tiny};

pub struct Quaternion<T> {
    pub c: [ T; 4 ]         // c[0] = a, c[1] = b, c[2] = c, c[3] = d 
}

pub type Quat64 = Quaternion<f64>;

impl<T> Quaternion<T> {
    /// Create a new quaternion ( q = a + bi + cj + dk )
    #[inline]
    pub const fn new(a: T, b: T, c: T, d: T) -> Self {
        Quaternion { c: [a, b, c, d] }
    }
}

impl<T: Clone + Copy + Signed> Quaternion<T> {
    /// Return the quaternion conjugate ( q.conj() = a - bi - cj - dk )
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(
            self.c[0], 
            -self.c[1], 
            -self.c[2], 
            -self.c[3]
        )
    }

    /// Return the inverse of the quaternion ( q^-1 = q.conj() / |q|^2 )
    #[inline]
    pub fn inv(&self) -> Self {
        let norm_sq = self.c[0] * self.c[0] + self.c[1] * self.c[1] 
                    + self.c[2] * self.c[2] + self.c[3] * self.c[3];
        Self::new(
            self.c[0] / norm_sq, 
            -self.c[1] / norm_sq, 
            -self.c[2] / norm_sq, 
            -self.c[3] / norm_sq
        )
    }

}

impl<T: Clone> Clone for Quaternion<T> {
    /// Clone the quaternion
    #[inline]
    fn clone(&self) -> Self {
        Self::new(
            self.c[0].clone(), 
            self.c[1].clone(), 
            self.c[2].clone(), 
            self.c[3].clone()
        )
    }
}

impl<T: Clone + Copy + Signed> Neg for Quaternion<T> {
    type Output = Self;
    /// Return the unary negation ( unary - )
    /// - ( a + bi + cj + dk ) = -a - bi - cj - dk
    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new( 
            -self.c[0], 
            -self.c[1], 
            -self.c[2], 
            -self.c[3] 
        )
    }
}

/* Addition */

// Full non-consuming addition
impl<T: Clone + Copy + Number> Add<&Quaternion<T>> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Add two quaternions together ( binary + )
    /// ( a + bi + cj + dk ) + ( e + fi + gj + hk ) 
    /// = ( a + e ) + i( b + f ) + j( c + g ) + k( d + h )
    #[inline]
    fn add(self, plus: &Quaternion<T>) -> Self::Output {
        Self::Output::new(
            self.c[0] + plus.c[0], 
            self.c[1] + plus.c[1], 
            self.c[2] + plus.c[2], 
            self.c[3] + plus.c[3]
        )
    }
}

// Partial non-consuming addition
impl<T: Clone + Copy + Number> Add<&Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Add two quaternions together ( binary + )
    #[inline]
    fn add(self, plus: &Self) -> Self::Output {
        &self + plus
    }
}

// Consuming addition
impl<T: Clone + Copy + Number> Add<Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Add two quaternions together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        self + &plus
    }
}

// Full non-consuming addition with a real number
impl<T: Clone + Copy + Number> Add<&T> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Add a real number to a quaternion
    /// ( a + bi + cj + dk ) + e = ( a + e ) + bi + cj + dk
    #[inline]
    fn add(self, plus: &T) -> Self::Output {
        Self::Output::new( self.c[0] + *plus, self.c[1], self.c[2], self.c[3] )
    }
}

// Partial non-consuming addition with a real number
impl<T: Clone + Copy + Number> Add<&T> for Quaternion<T> {
    type Output = Self;
    /// Add a real number to a quaternion
    #[inline]    
    fn add(self, plus: &T) -> Self::Output {
        &self + plus
    }
}

// Consuming addition with a real number
impl<T: Clone + Copy + Number> Add<T> for Quaternion<T> {
    type Output = Self;
    /// Add a real number to a quaternion
    #[inline]    
    fn add(self, plus: T) -> Self::Output {
        self + &plus
    }
}

/* Subtraction */

// Full non-consuming subtraction
impl<T: Clone + Copy + Number> Sub<&Quaternion<T>> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Subtract one quaternion from another ( binary - )
    /// ( a + bi + cj + dk ) - ( e + fi + gj + hk ) 
    /// = ( a - e ) + i( b - f ) + j( c - g ) + k( d - h )
    #[inline]
    fn sub(self, minus: &Quaternion<T>) -> Self::Output {
        Self::Output::new(
            self.c[0] - minus.c[0], 
            self.c[1] - minus.c[1], 
            self.c[2] - minus.c[2], 
            self.c[3] - minus.c[3]
        )
    }
}

// Partial non-consuming subtraction
impl<T: Clone + Copy + Number> Sub<&Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Subtract one quaternion from another ( binary - )
    #[inline]
    fn sub(self, minus: &Self) -> Self::Output {
        &self - minus
    }
}

// Consuming subtraction
impl<T: Clone + Copy + Number> Sub<Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Subtract one quaternion from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        self - &minus
    }
}

// Full non-consuming subtraction with a real number
impl<T: Clone + Copy + Number> Sub<&T> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Subtract a real number from a quaternion
    /// ( a + bi + cj + dk ) - e = ( a - e ) + bi + cj + dk
    #[inline]
    fn sub(self, minus: &T) -> Self::Output {
        Self::Output::new( self.c[0] - *minus, self.c[1], self.c[2], self.c[3] )
    }
}

// Partial non-consuming subtraction with a real number
impl<T: Clone + Copy + Number> Sub<&T> for Quaternion<T> {
    type Output = Self;
    /// Subtract a real number from a quaternion
    #[inline]    
    fn sub(self, minus: &T) -> Self::Output {
        &self - minus
    }
}

// Consuming subtraction with a real number
impl<T: Clone + Copy + Number> Sub<T> for Quaternion<T> {
    type Output = Self;
    /// Subtract a real number from a quaternion
    #[inline]    
    fn sub(self, minus: T) -> Self::Output {
        self - &minus
    }
}

/* Multiplication */

// Full non-consuming multiplication
impl<T: Clone + Copy + Number> Mul<&Quaternion<T>> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Multiply two quaternions together ( binary * )
    /// ( a + bi + cj + dk ) * ( e + fi + gj + hk ) 
    /// = ( ae - bf - cg - dh ) + i( af + be + ch - dg ) 
    /// + j( ag - bh + ce + df ) + k( ah + bg - cf + de )
    #[inline]
    fn mul(self, times: &Quaternion<T>) -> Self::Output {
        Self::Output::new(
            self.c[0] * times.c[0] - self.c[1] * times.c[1] - self.c[2] * times.c[2] - self.c[3] * times.c[3],
            self.c[0] * times.c[1] + self.c[1] * times.c[0] + self.c[2] * times.c[3] - self.c[3] * times.c[2],
            self.c[0] * times.c[2] - self.c[1] * times.c[3] + self.c[2] * times.c[0] + self.c[3] * times.c[1],
            self.c[0] * times.c[3] + self.c[1] * times.c[2] - self.c[2] * times.c[1] + self.c[3] * times.c[0]
        )
    }
}

// Partial non-consuming multiplication
impl<T: Clone + Copy + Number> Mul<&Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Multiply two quaternions together ( binary * )
    #[inline]
    fn mul(self, times: &Self) -> Self::Output {
        &self * times
    }
}


// Consuming multiplication
impl<T: Clone + Copy + Number> Mul<Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Multiply two quaternions together ( binary * )
    #[inline]
    fn mul(self, times: Self) -> Self::Output {
        self * &times
    }
}

// Full non-consuming multiplication by a real number
impl<T: Clone + Copy + Number> Mul<&T> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Multiply a quaternion by a real number
    /// ( a + bi + cj + dk ) * e = ae + bei + cej + dek
    #[inline]
    fn mul(self, times: &T) -> Self::Output {
        Self::Output::new( 
            self.c[0] * *times, 
            self.c[1] * *times, 
            self.c[2] * *times, 
            self.c[3] * *times
        )
    }
}

// Partial non-consuming multiplication by a real number
impl<T: Clone + Copy + Number> Mul<&T> for Quaternion<T> {
    type Output = Self;
    /// Multiply a quaternion by a real number
    #[inline]    
    fn mul(self, times: &T) -> Self::Output {
        &self * times
    }
}

// Consuming multiplication by a real number
impl<T: Clone + Copy + Number> Mul<T> for Quaternion<T> {
    type Output = Self;
    /// Multiply a quaternion by a real number
    #[inline]    
    fn mul(self, times: T) -> Self::Output {
        self * &times
    }
}

impl Mul<Quaternion<f64>> for f64 {
    type Output = Quaternion<f64>;
    /// Allow multiplication on the left by f64 (f64 * quaternion number)
    #[inline]
    fn mul(self, q: Quaternion<f64>) -> Self::Output {
        q * self
    }
}

/* Division */

// Full non-consuming division
impl<T: Clone + Copy + Number> Div<&Quaternion<T>> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Divide one quaternion by another ( binary / )
    /// ( a + bi + cj + dk ) / ( e + fi + gj + hk ) 
    /// = [( ae + bf + cg + dh ) + i( - af + be - ch + dg ) 
    /// + j( - ag + bh + ce - df ) + k( - ah - bg + cf + de )] / ( e^2 + f^2 + g^2 + h^2 )
    #[inline]
    fn div(self, divisor: &Quaternion<T>) -> Self::Output {
        let norm_sq = divisor.c[0] * divisor.c[0] + divisor.c[1] * divisor.c[1] 
                    + divisor.c[2] * divisor.c[2] + divisor.c[3] * divisor.c[3];
        Self::Output::new(
            (self.c[0] * divisor.c[0] + self.c[1] * divisor.c[1] 
                + self.c[2] * divisor.c[2] + self.c[3] * divisor.c[3]) / norm_sq,
            (self.c[1] * divisor.c[0] - self.c[0] * divisor.c[1] 
                - self.c[2] * divisor.c[3] + self.c[3] * divisor.c[2]) / norm_sq,
            (self.c[2] * divisor.c[0] - self.c[0] * divisor.c[2] 
                + self.c[1] * divisor.c[3] - self.c[3] * divisor.c[1]) / norm_sq,
            (self.c[3] * divisor.c[0] - self.c[0] * divisor.c[3] 
                - self.c[1] * divisor.c[2] + self.c[2] * divisor.c[1]) / norm_sq
        )
    }
}

// Partial non-consuming division
impl<T: Clone + Copy + Number> Div<&Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Divide one quaternion by another ( binary / )
    #[inline]    
    fn div(self, divisor: &Self) -> Self::Output {
        &self / divisor
    }
}

// Consuming division
impl<T: Clone + Copy + Number> Div<Quaternion<T>> for Quaternion<T> {
    type Output = Self;
    /// Divide one quaternion by another ( binary / )
    #[inline]    
    fn div(self, divisor: Self) -> Self::Output {
        self / &divisor
    }
}

// Full non-consuming division by a real number
impl<T: Clone + Copy + Number> Div<&T> for &Quaternion<T> {
    type Output = Quaternion<T>;
    /// Divide a quaternion by a real number
    /// ( a + bi + cj + dk ) / e = (a/e) + (b/e)i + (c/e)j + (d/e)k
    #[inline]
    fn div(self, divisor: &T) -> Self::Output {
        Self::Output::new( 
            self.c[0] / *divisor, 
            self.c[1] / *divisor, 
            self.c[2] / *divisor, 
            self.c[3] / *divisor 
        )
    }
}

// Partial non-consuming division by a real number
impl<T: Clone + Copy + Number> Div<&T> for Quaternion<T> {
    type Output = Self;
    /// Divide a quaternion by a real number
    #[inline]    
    fn div(self, divisor: &T) -> Self::Output {
        &self / divisor
    }
}

// Consuming division by a real number
impl<T: Clone + Copy + Number> Div<T> for Quaternion<T> {
    type Output = Self;
    /// Divide a quaternion by a real number
    #[inline]    
    fn div(self, divisor: T) -> Self::Output {
        self / &divisor
    }
}

/* Add assign */

// Non-consuming add assign
impl<T: Clone + Copy + Number> AddAssign<&Quaternion<T>> for Quaternion<T> {
    /// Add-assign one quaternion to another ( binary += )
    /// ( a + bi + cj + dk ) += ( e + fi + gj + hk ) 
    /// = ( a + e ) + i( b + f ) + j( c + g ) + k( d + h )
    #[inline]
    fn add_assign(&mut self, plus: &Quaternion<T>) {
        self.c[0] += plus.c[0];
        self.c[1] += plus.c[1];
        self.c[2] += plus.c[2];
        self.c[3] += plus.c[3];
    }
}

// Consuming add assign
impl<T: Clone + Copy + Number> AddAssign<Quaternion<T>> for Quaternion<T> {
    /// Add-assign one quaternion to another ( binary += )
    #[inline]
    fn add_assign(&mut self, plus: Quaternion<T>) {
        *self += &plus;
    }
}

// Non-consuming add assign with a real number
impl<T: Clone + Copy + Number> AddAssign<&T> for Quaternion<T> {
    /// Add-assign a real number to a quaternion
    /// ( a + bi + cj + dk ) += e = ( a + e ) + bi + cj + dk
    #[inline]
    fn add_assign(&mut self, plus: &T) {
        self.c[0] += *plus;
    }
}

// Consuming add assign with a real number
impl<T: Clone + Copy + Number> AddAssign<T> for Quaternion<T> {
    /// Add-assign a real number to a quaternion
    #[inline]
    fn add_assign(&mut self, plus: T) {
        *self += &plus;
    }
}

/* Subtract assign */

// Non-consuming sub assign
impl<T: Clone + Copy + Number> SubAssign<&Quaternion<T>> for Quaternion<T> {
    /// Subtract-assign one quaternion from another ( binary -= )
    /// ( a + bi + cj + dk ) -= ( e + fi + gj + hk ) 
    /// = ( a - e ) + i( b - f ) + j( c - g ) + k( d - h )
    #[inline]
    fn sub_assign(&mut self, minus: &Quaternion<T>) {
        self.c[0] -= minus.c[0];
        self.c[1] -= minus.c[1];
        self.c[2] -= minus.c[2];
        self.c[3] -= minus.c[3];
    }
}

// Consuming sub assign
impl<T: Clone + Copy + Number> SubAssign<Quaternion<T>> for Quaternion<T> {
    /// Subtract-assign one quaternion from another ( binary -= )
    #[inline]
    fn sub_assign(&mut self, minus: Quaternion<T>) {
        *self -= &minus;
    }
}

// Non-consuming sub assign with a real number
impl<T: Clone + Copy + Number> SubAssign<&T> for Quaternion<T> {
    /// Subtract-assign a real number from a quaternion
    /// ( a + bi + cj + dk ) -= e = ( a - e ) + bi + cj + dk
    #[inline]
    fn sub_assign(&mut self, minus: &T) {
        self.c[0] -= *minus;
    }
}

// Consuming sub assign with a real number
impl<T: Clone + Copy + Number> SubAssign<T> for Quaternion<T> {
    /// Subtract-assign a real number from a quaternion
    #[inline]
    fn sub_assign(&mut self, minus: T) {
        *self -= &minus;
    }
}

/* Multiply assign */

// Non-consuming mul assign
impl<T: Clone + Copy + Number> MulAssign<&Quaternion<T>> for Quaternion<T> {
    /// Multiply-assign one quaternion by another ( binary *= )
    /// ( a + bi + cj + dk ) *= ( e + fi + gj + hk ) 
    /// = ( ae - bf - cg - dh ) + i( af + be + ch - dg ) 
    /// + j( ag - bh + ce + df ) + k( ah + bg - cf + de )
    #[inline]
    fn mul_assign(&mut self, times: &Quaternion<T>) {
        let [a, b, c, d] = self.c;
        self.c[0] = a * times.c[0] - b * times.c[1] - c * times.c[2] - d * times.c[3];
        self.c[1] = a * times.c[1] + b * times.c[0] + c * times.c[3] - d * times.c[2];
        self.c[2] = a * times.c[2] - b * times.c[3] + c * times.c[0] + d * times.c[1];
        self.c[3] = a * times.c[3] + b * times.c[2] - c * times.c[1] + d * times.c[0];
    }
}

// Consuming mul assign
impl<T: Clone + Copy + Number> MulAssign<Quaternion<T>> for Quaternion<T> {
    /// Multiply-assign one quaternion by another ( binary *= )
    #[inline]
    fn mul_assign(&mut self, times: Quaternion<T>) {
        *self *= &times;
    }
}

// Non-consuming mul assign by a real number
impl<T: Clone + Copy + Number> MulAssign<&T> for Quaternion<T> {
    /// Multiply-assign a quaternion by a real number
    /// ( a + bi + cj + dk ) *= e = ae + bei + cej + dek
    #[inline]
    fn mul_assign(&mut self, times: &T) {
        self.c[0] *= *times;
        self.c[1] *= *times;
        self.c[2] *= *times;
        self.c[3] *= *times;
    }
}

// Consuming mul assign by a real number
impl<T: Clone + Copy + Number> MulAssign<T> for Quaternion<T> {
    /// Multiply-assign a quaternion by a real number
    #[inline]    
    fn mul_assign(&mut self, times: T) {
        *self *= &times;
    }
}

/* Divide assign */

// Non-consuming div assign
impl<T: Clone + Copy + Number> DivAssign<&Quaternion<T>> for Quaternion<T> {
    /// Divide-assign one quaternion by another ( binary /= )
    /// ( a + bi + cj + dk ) /= ( e + fi + gj + hk ) 
    /// = [( ae + bf + cg + dh ) + i( - af + be - ch + dg ) 
    /// + j( - ag + bh + ce - df ) + k( - ah - bg + cf + de )] / ( e^2 + f^2 + g^2 + h^2 )
    #[inline]
    fn div_assign(&mut self, div: &Quaternion<T>) {
        let norm_sq = div.c[0] * div.c[0] + div.c[1] * div.c[1] 
                    + div.c[2] * div.c[2] + div.c[3] * div.c[3];
        let [a, b, c, d] = self.c;
        self.c[0] = (a * div.c[0] + b * div.c[1] + c * div.c[2] + d * div.c[3]) / norm_sq;
        self.c[1] = (b * div.c[0] - a * div.c[1] - c * div.c[3] + d * div.c[2]) / norm_sq;
        self.c[2] = (c * div.c[0] - a * div.c[2] + b * div.c[3] - d * div.c[1]) / norm_sq;
        self.c[3] = (d * div.c[0] - a * div.c[3] - b * div.c[2] + c * div.c[1]) / norm_sq;
    }
}

// Consuming div assign
impl<T: Clone + Copy + Number> DivAssign<Quaternion<T>> for Quaternion<T> {
    /// Divide-assign one quaternion by another ( binary /= )
    #[inline]
    fn div_assign(&mut self, div: Quaternion<T>) {
        *self /= &div;
    }
}

// Non-consuming div assign by a real number
impl<T: Clone + Copy + Number> DivAssign<&T> for Quaternion<T> {
    /// Divide-assign a quaternion by a real number
    /// ( a + bi + cj + dk ) /= e = (a/e) + (b/e)i + (c/e)j + (d/e)k
    #[inline]
    fn div_assign(&mut self, div: &T) {
        self.c[0] /= *div;
        self.c[1] /= *div;
        self.c[2] /= *div;
        self.c[3] /= *div;
    }
}

// Consuming div assign by a real number
impl<T: Clone + Copy + Number> DivAssign<T> for Quaternion<T> {
    /// Divide-assign a quaternion by a real number
    #[inline]    
    fn div_assign(&mut self, div: T) {
        *self /= &div;
    }
}

/* Numbers */

impl<T: Clone + Copy + Number> Zero for Quaternion<T> {
    /// Return the additive identity q = 0 + 0i + 0j + 0k 
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero(), Zero::zero(), Zero::zero())
    }
}

impl<T: Clone + Copy + Number> One for Quaternion<T> {
    /// Return the multiplicative identity q = 1 + 0i + 0j + 0k
    fn one() -> Self {
        Self::new(One::one(), Zero::zero(), Zero::zero(), Zero::zero())
    }
}

impl<T: Clone + Copy + Number> Tiny for Quaternion<T> {
    /// Return a very small quaternion number q = tiny + tiny*i + tiny*j + tiny*k
    fn tiny() -> Self {
        Self::new(T::tiny(), T::tiny(), T::tiny(), T::tiny())
    }
}

/* Comparison */

impl<T: Clone + Number> PartialEq for Quaternion<T> {
    /// Implement trait for equality 
    fn eq(&self, q: &Quaternion<T>) -> bool {
        self.c[0] == q.c[0] && self.c[1] == q.c[1] && self.c[2] == q.c[2] && self.c[3] == q.c[3]
    }
}

impl<T: Clone + Number + std::cmp::PartialOrd> PartialOrd for Quaternion<T> {
    /// Implement trait for ordering 
    fn partial_cmp(&self, other: &Quaternion<T>) -> Option<Ordering> {
        if self.c[0] != other.c[0] {
            self.c[0].partial_cmp( &other.c[0] )
        } else {
            self.c[1..].partial_cmp( &other.c[1..] )
        }
    }
}

/* Functions */

impl<T: Clone + Copy + Number> Quaternion<T> {
    /// Return the absolute value squared ( |q|^2 = q * q.conj() )
    #[inline]
    pub fn abs_sqr(&self) -> T {
        self.c[0] * self.c[0] + self.c[1] * self.c[1] + 
        self.c[2] * self.c[2] + self.c[3] * self.c[3]
    }

    /// Return the scalar part of the quaternion ( a in q = a + bi + cj + dk )
    #[inline]    
    pub fn scalar(&self) -> T {
        self.c[0]
    }

    /// Return the vector part of the quaternion ( v = bi + cj + dk )
    #[inline]
    pub fn vector(&self) -> [T; 3] {
        [self.c[1], self.c[2], self.c[3]]
    }
}

impl Quaternion<f64> {
    /// Return the norm ( |q| = sqrt( q * q.conj() ) )
    #[inline]
    pub fn norm(&self) -> f64 {
        f64::sqrt( self.abs_sqr() )
    }

    /// Return the versor ( unit quaternion ( q / |q| ) )
    #[inline]
    pub fn versor(&self) -> Self {
        let norm = self.norm();
        Self::new(
            self.c[0] / norm, 
            self.c[1] / norm, 
            self.c[2] / norm, 
            self.c[3] / norm
        )
    }
}

impl<T> fmt::Display for Quaternion<T> where
    T: fmt::Display
{
    /// Format the output ( z = a + bi + cj + dk )
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}i + {}j + {}k", self.c[0], self.c[1], self.c[2], self.c[3])
    }
}

impl<T> fmt::Debug for Quaternion<T> where
    T: fmt::Debug
{
    /// Format the output ( z = a + bi + cj + dk )
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Quaternion {{ c: {:?} }}", self.c )
    }
} 

impl<T: Number + Copy + Clone> Number for Quaternion<T> {

}

impl Signed for Quaternion<f64> {
    fn abs(&self) -> Self {
        let a: f64 = self.norm();
        Quaternion { c: [a, 0.0, 0.0, 0.0] }
    }
}

impl Copy for Quaternion<f64> {

}