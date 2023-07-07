use crate::complex::{Complex, Cmplx};
use crate::traits::One;

impl Complex::<f64> {
    /// Return the hyperbolic sine of a complex number z ( sinh(z) )
    #[inline]
    pub fn sinh(&self) -> Complex::<f64> {
        Cmplx::new(self.real.sinh() * self.imag.cos(), self.real.cosh() * self.imag.sin())
    }

    /// Return the hyperbolic cosine of a complex number z ( cosh(z) )
    #[inline]
    pub fn cosh(&self) -> Complex::<f64> {
        Cmplx::new(self.real.cosh() * self.imag.cos(), self.real.sinh() * self.imag.sin())
    }

    /// Return the hyperbolic tangent of a complex number z ( tanh(z) )
    #[inline]
    pub fn tanh(&self) -> Complex::<f64> {
        self.sinh() / self.cosh()
    }

    /// Return the hyperbolic secant of a complex number z ( sech(z) )
    #[inline]
    pub fn sech(&self) -> Complex::<f64> {
        Cmplx::one() / self.cosh()
    }

    /// Return the hyperbolic cosecant of a complex number z ( csch(z) )
    #[inline]
    pub fn csch(&self) -> Complex::<f64> {
        Cmplx::one() / self.sinh()
    }

    /// Return the hyperbolic cotangent of a complex number z ( coth(z) )
    #[inline]
    pub fn coth(&self) -> Complex::<f64> {
        Cmplx::one() / self.tanh()
    }

    /// Return the inverse hyperbolic sine of a complex number z ( asinh(z) )
    #[inline]
    pub fn asinh(&self) -> Complex::<f64> {
        let z = self.clone();
        ( ( z * z + 1.0 ).sqrt() + z ).ln()
    }

    /// Return the inverse hyperbolic cosine of a complex number z ( acosh(z) )
    #[inline]
    pub fn acosh(&self) -> Complex::<f64> {
        let z = self.clone();
        ( ( z - 1.0 ).sqrt() * ( z + 1.0 ).sqrt() + z ).ln()
    }

    /// Return the inverse hyperbolic tangent of a complex number z ( atanh(z) )
    #[inline]
    pub fn atanh(&self) -> Complex::<f64> {
        let z = self.clone();
        ( ( z + 1.0 ).ln() - ( Cmplx::one() - z ).ln() ) * 0.5
    }

    /// Return the inverse hyperbolic secant of a complex number z ( asech(z) )
    #[inline]
    pub fn asech(&self) -> Complex::<f64> {
        let inv = Cmplx::one() / self.clone();
        inv.acosh()
    }

    /// Return the inverse hyperbolic cosecant of a complex number z ( acsch(z) )
    #[inline]
    pub fn acsch(&self) -> Complex::<f64> {
        let inv = Cmplx::one() / self.clone();
        inv.asinh()
    }
    
    /// Return the inverse hyperbolic cotangent of a complex number z ( acoth(z) )
    #[inline]
    pub fn acoth(&self) -> Complex::<f64> {
        let inv = Cmplx::one() / self.clone();
        inv.atanh()
    }
}
