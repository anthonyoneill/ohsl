use core::ops::{Add, Div, Mul, Neg, Sub};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::{fmt, cmp::Ordering};

pub use crate::traits::{Number, Signed, Zero, One};

pub struct Complex<T> {
    /// Real part of the complex number
    pub real: T,
    /// Imaginary part of the complex number
    pub imag: T,
}

pub type Cmplx = Complex<f64>;

impl<T> Complex<T> {
    /// Create a new complex number ( z = x + iy )
    #[inline]
    pub const fn new(real: T, imag: T) -> Self {
        Complex { real, imag }
    }
}

impl<T: Clone + Signed> Complex<T> {
    /// Return the complex conjugate ( z.conj() = x - iy )
    #[inline]
    pub fn conj(&self) -> Self {
        Self::new(self.real.clone(), -self.imag.clone())
    }
}

impl<T: Clone> Clone for Complex<T> {
    /// Clone the complex number
    #[inline]
    fn clone(&self) -> Self {
        Self::new(self.real.clone(), self.imag.clone())
    }
}

impl<T: Clone + Signed> Neg for Complex<T> {
    type Output = Self;
    /// Return the unary negation ( unary - )
    /// - ( a + ib ) = -a - ib
    #[inline]
    fn neg(self) -> Self::Output {
        Self::Output::new( -self.real, -self.imag )
    }
}

impl<T: Clone + Number> Add<Complex<T>> for Complex<T> {
    type Output = Self;
    /// Add two complex numbers together ( binary + )
    /// ( a + ib ) + ( c + id ) = ( a + c ) + i( b + d )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        Self::Output::new(self.real + plus.real, self.imag + plus.imag)
    }
}

impl<T: Clone + Number> Sub<Complex<T>> for Complex<T> {
    type Output = Self;
    /// Subtract one complex number from another ( binary - )
    /// ( a + ib ) - ( c + id ) = ( a - c ) + i( b - d )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        Self::Output::new(self.real - minus.real, self.imag - minus.imag)
    }
}

impl<T: Clone + Number> Mul<Complex<T>> for Complex<T> {
    type Output = Self;
    /// Multiply two complex numbers together ( binary * )
    /// ( a + ib ) * ( c + id ) = ( ac - bd ) + i( ad + bc )
    #[inline]
    fn mul(self, times: Self) -> Self::Output {
        let real = self.real.clone() * times.real.clone() - self.imag.clone() * times.imag.clone();  
        let imag = self.real * times.imag + self.imag * times.real;
        Self::Output::new( real, imag )
    }
}

impl<T: Clone + Number> Div<Complex<T>> for Complex<T> {
    type Output = Self;
    /// Divide one complex number by another ( binary / )
    /// ( a + ib ) / ( c + id ) = [( ac + bd ) + i( bc - ad )] / ( c^2 + d^2 )
    #[inline]
    fn div(self, divisor: Self) -> Self::Output {
        let denominator = divisor.real.clone() * divisor.real.clone() + divisor.imag.clone() * divisor.imag.clone();
        let real = self.real.clone() * divisor.real.clone() + self.imag.clone() * divisor.imag.clone();  
        let imag = self.imag * divisor.real - self.real * divisor.imag;
        Self::Output::new( real / denominator.clone(), imag / denominator )
    }
}

impl<T: Number> Add<T> for Complex<T> {
    type Output = Complex<T>;
    /// Add a complex number to a real number 
    /// ( a + ib ) + c = ( a + c ) + ib 
    fn add(self, plus: T) -> Self::Output {
        Self::Output::new( self.real + plus, self.imag )
    }
}

impl<T: Number> Sub<T> for Complex<T> {
    type Output = Complex<T>;
    /// Subtract a real number from a complex number 
    /// ( a + ib ) - c = ( a - c ) + ib 
    fn sub(self, minus: T) -> Self::Output {
        Self::Output::new( self.real - minus, self.imag )
    }
}

impl<T: Clone + Number> Mul<T> for Complex<T>  {
    type Output = Complex<T>;
    /// Multiply a complex number by a real scalar 
    /// ( a + ib ) * r = (a*r) + i(b*r) 
    fn mul(self, scalar: T) -> Self::Output {
        Self::Output::new( self.real * scalar.clone(), self.imag * scalar )
    }
}

impl<T: Clone + Number> Div<T> for Complex<T> {
    type Output = Complex<T>;
    /// Divide a complex number by a real scalar 
    /// ( a + ib ) / r = (a/r) + i(b/r) 
    fn div(self, scalar: T) -> Self::Output {
        Self::Output::new( self.real / scalar.clone(), self.imag / scalar )
    }
}

impl<T: Number> AddAssign for Complex<T> {
    /// Add a complex number to a mutable complex variable and assign the
    /// result to that variable ( += )
    fn add_assign(&mut self, rhs: Self) {
        self.real += rhs.real;
        self.imag += rhs.imag;
    }
}

impl<T: Number> SubAssign for Complex<T> {
    /// Subtract a complex number from a mutable complex variable and assign the
    /// result to that variable ( -= )
    fn sub_assign(&mut self, rhs: Self) {
        self.real -= rhs.real;
        self.imag -= rhs.imag;
    }
}

impl<T: Clone + Number> MulAssign for Complex<T> {
    /// Multipy a mutable complex variable by a complex number and assign the
    /// result to that variable ( *= )
    fn mul_assign(&mut self, rhs: Self) {
        let a = self.real.clone();

        self.real *= rhs.real.clone();
        self.real -= self.imag.clone() * rhs.imag.clone();

        self.imag *= rhs.real;
        self.imag += a * rhs.imag;
    }
}

impl<T: Clone + Number> DivAssign for Complex<T> {
    /// Divide a mutable complex variable by a complex number and assign the
    /// result to that variable ( *= )
    fn div_assign(&mut self, rhs: Self) {
        let a = self.real.clone();
        let denominator = rhs.real.clone() * rhs.real.clone() + rhs.imag.clone() * rhs.imag.clone();
        
        self.real *= rhs.real.clone();
        self.real += self.imag.clone() * rhs.imag.clone();
        self.real /= denominator.clone();

        self.imag *= rhs.real;
        self.imag -= a * rhs.imag;
        self.imag /= denominator;
    }
}

impl<T: Number> AddAssign<T> for Complex<T> {
    /// Add a real number to a mutable complex variable and assign the
    /// result to that variable ( += )
    fn add_assign(&mut self, rhs: T) {
        self.real += rhs;
    }
}

impl<T: Number> SubAssign<T> for Complex<T> {
    /// Subtract a real number from a mutable complex variable and assign the
    /// result to that variable ( += )
    fn sub_assign(&mut self, rhs: T) {
        self.real -= rhs;
    }
}

impl<T: Clone + Number> MulAssign<T> for Complex<T> {
    /// Multiply a mutable complex variable by a real number and assign the
    /// result to that variable ( *= )
    fn mul_assign(&mut self, rhs: T) {
        self.real *= rhs.clone();
        self.imag *= rhs;
    }
}

impl<T: Clone + Number> DivAssign<T> for Complex<T> {
    /// Divide a mutable complex variable by a real number and assign the
    /// result to that variable ( /= )
    fn div_assign(&mut self, rhs: T) {
        self.real /= rhs.clone();
        self.imag /= rhs;
    }
}

impl<T: Clone + Number> Zero for Complex<T> {
    /// Return the additive identity z = 0 + 0i
    fn zero() -> Self {
        Self::new(Zero::zero(), Zero::zero())
    }
}

impl<T: Clone + Number> One for Complex<T> {
    /// Return the multiplicative identity z = 1 + 0i
    fn one() -> Self {
        Self::new(One::one(), Zero::zero())
    }
}

impl<T: Clone + Number> PartialEq for Complex<T> {
    /// Implement trait for equality 
    fn eq(&self, other: &Complex<T>) -> bool {
        self.real == other.real && self.imag == other.imag
    }
}

impl<T: Clone + Number + std::cmp::PartialOrd> PartialOrd for Complex<T> {
    /// Implement trait for ordering 
    fn partial_cmp(&self, other: &Complex<T>) -> Option<Ordering> {
        if self.real != other.real {
            self.real.partial_cmp( &other.real )
        } else {
            self.imag.partial_cmp( &other.imag )
        }
    }
}

impl<T: Clone + Number> Complex<T> {
    /// Return the absolute value squared ( |z|^2 = z * z.conj() )
    #[inline]
    pub fn abs_sqr(&self) -> T {
        self.real.clone() * self.real.clone() + self.imag.clone() * self.imag.clone()
    }
}

impl Complex::<f64> {
    /// Return the absolute value ( |z| = sqrt( z * z.conj() ) )
    #[inline]
    pub fn abs(&self) -> f64 {
        libm::sqrt( self.abs_sqr() )
    }
}

impl Complex::<f64> {
    /// Return the phase angle in radians
    #[inline]
    pub fn arg(&self) -> f64 {
        self.imag.atan2(self.real)
    }
}

impl<T> fmt::Display for Complex<T> where
    T: fmt::Display
{
    /// Format the output ( z = x + iy = ( x, y ) )
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "( {}, {} )", self.real, self.imag )
    }
} 

impl<T> fmt::Debug for Complex<T> where
    T: fmt::Debug
{
    /// Format the output ( z = x + iy = ( x, y ) )
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "( {:?}, {:?} )", self.real, self.imag )
    }
} 

impl<T: Number + Clone> Number for Complex<T> {

}

impl Signed for Complex<f64> {
    fn abs(&self) -> Self {
        //let real = self.abs();
        let real: f64 = self.abs();
        let imag: f64 = 0.0;
        Complex { real, imag }
    }
}

impl Copy for Complex<f64> {

}