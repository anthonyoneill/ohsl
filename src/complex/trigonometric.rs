use crate::complex::{Complex, Cmplx};
use crate::traits::One;
use crate::constant::{I, PI_2};

impl Complex::<f64> {
    ///Return the sin of a complex number z ( sin(z) )
    #[inline]
    pub fn sin(&self) -> Complex<f64> {
        Cmplx::new(self.real.sin() * self.imag.cosh(), self.real.cos() * self.imag.sinh())
    }

    ///Return the cos of a complex number z ( cos(z) )
    #[inline]
    pub fn cos(&self) -> Complex<f64> {
        Cmplx::new(self.real.cos() * self.imag.cosh(), -self.real.sin() * self.imag.sinh())
    }

    ///Return the tan of a complex number z ( tan(z) )
    #[inline]
    pub fn tan(&self) -> Complex<f64> {
        self.sin() / self.cos()
    }

    ///Return the sec of a complex number z ( sec(z) )
    #[inline]
    pub fn sec(&self) -> Complex<f64> {
        Complex::<f64>::one() / self.cos()
    }

    ///Return the csc of a complex number z ( csc(z) )
    #[inline]
    pub fn csc(&self) -> Complex<f64> {
        Complex::<f64>::one() / self.sin()
    }
    
    ///Return the cot of a complex number z ( cot(z) )
    #[inline]
    pub fn cot(&self) -> Complex<f64> {
        Complex::<f64>::one() / self.tan()
    }

    /// Return the inverse sin of a complex number z ( asin(z) )
    #[inline]
    pub fn asin(&self) -> Complex<f64> {
        let squared = self.clone() * self.clone();
        - I * ((Cmplx::one() - squared).sqrt() + I * self.clone()).ln()
    }

    /// Return the inverse cos of a complex number z ( acos(z) )
    #[inline]
    pub fn acos(&self) -> Complex<f64> {
        let squared = self.clone() * self.clone();
        I * ((Cmplx::one() - squared).sqrt() + I * self.clone()).ln() + PI_2
    }

    /// Return the inverse tan of a complex number z ( atan(z) )
    #[inline]
    pub fn atan(&self) -> Complex<f64> {
        let iz = I * self.clone();
        ( (Cmplx::one() - iz).ln() - (Cmplx::one() + iz).ln() ) * I * 0.5
    }

    /// Return the inverse sec of a complex number z ( asec(z) )
    #[inline]
    pub fn asec(&self) -> Complex<f64> {
        let inv = Complex::<f64>::one() / self.clone();
        inv.acos()
    }

    /// Return the inverse csc of a complex number z ( acsc(z) )
    #[inline]
    pub fn acsc(&self) -> Complex<f64> {
        let inv = Complex::<f64>::one() / self.clone();
        inv.asin()
    }

    /// Return the inverse cot of a complex number z ( acot(z) )
    #[inline]
    pub fn acot(&self) -> Complex<f64> {
        let inv = Complex::<f64>::one() / self.clone();
        inv.atan()
    }
}
