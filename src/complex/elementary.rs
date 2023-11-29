use crate::complex::Complex;

impl Complex::<f64> {
    /// Return the square root of a complex number ( sqrt(z) )
    #[inline]
    pub fn sqrt(&self) -> Complex::<f64> {
        let sqrt_abs = (self.abs()).sqrt();
        let theta = self.arg();
        let x = sqrt_abs * f64::cos( 0.5 * theta );
        let y = sqrt_abs * f64::sin( 0.5 * theta );
        Complex::new( x, y )
    }

    /// Return the complex number z raised to the power of a complex number w ( z^w )
    #[inline]
    pub fn pow(&self, w: &Complex::<f64>) -> Complex::<f64> {
        let r2 = self.abs_sqr();
        let theta = self.arg();
        let x = r2.powf( 0.5 * w.real ) * f64::exp( -w.imag * theta );
        let y = w.real * theta + 0.5 * w.imag * f64::ln( r2 );
        Complex::new( x * f64::cos(y), x * f64::sin(y) )
    }

    /// Return complex number z raised to the power of a real number x ( z^x )
    #[inline]
    pub fn powf(&self, x: f64) -> Complex::<f64> {
        let r2 = self.abs_sqr();
        let theta = self.arg();
        let a = r2.powf( 0.5 * x );
        let b = x * theta;
        Complex::new( a * f64::cos(b), a * f64::sin(b) )
    }

    /// Return the complex exponential of the complex number z ( exp(z) )
    #[inline]
    pub fn exp(&self) -> Complex::<f64> {
        let a = f64::exp( self.real );
        Complex::new( a * f64::cos(self.imag), a * f64::sin(self.imag) )
    }

    /// Return the complex logarithm of the complex number z ( ln(z) )
    #[inline]
    pub fn ln(&self) -> Complex::<f64> {
        let r = self.abs();
        let theta = self.arg();
        Complex::new( f64::ln(r), theta )
    }

    /// Return the complex base-b logarithm of the complex number z ( log_b(z) )
    #[inline]
    pub fn log(&self, b: Complex::<f64>) -> Complex::<f64> {
        self.ln() / b.ln()
    }

    /// Create a new complex number with magnitude r and phase angle theta
    #[inline]
    pub fn polar(r: f64, theta: f64) -> Complex::<f64> {
        let real = r * theta.cos();
        let imag = r * theta.sin();
        Complex::new( real, imag )
    }
}