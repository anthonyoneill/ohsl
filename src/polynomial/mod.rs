pub mod arithmetic;

use core::ops::{Add, Mul};
use std::fmt;
use crate::{ Cmplx, Vector };
use crate::traits::{Zero, Signed};

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

    /// Create a new quadratic polynomial of the form ax^2 + bx + c
    #[inline]
    pub fn quadratic( a: T, b: T, c: T ) -> Self {
        Polynomial { coeffs: vec![ c, b, a ] }
    }

    /// Create a new cubic polynomial of the form ax^3 + bx^2 + cx + d
    #[inline]
    pub fn cubic( a: T, b: T, c: T, d: T ) -> Self {
        Polynomial { coeffs: vec![ d, c, b, a ] }
    }
    
    /// Return the size of the polynomial
    #[inline]
    pub fn size(&self) -> usize {
        self.coeffs.len()
    }

    /// Return the degree of the polynomial or an error if < 1.
    #[inline]
    pub fn degree(&self) -> Result<usize, &'static str> {
        if self.coeffs.len() == 0 { Err("Polynomial.degree() == 0.") }
        else { Ok( self.coeffs.len() - 1 ) }
    }

    /// Return the coefficients of the polynomial as a mutable reference
    #[inline]
    pub fn coeffs(&mut self) -> &mut Vec<T> {
        &mut self.coeffs
    }

    /// Evaluate the polynomial at a given point
    #[inline]
    pub fn eval(&self, x: T) -> T
    where
        T: Copy + Mul<Output = T> + Add<Output = T>,
    {
        let degree = self.degree().unwrap(); //TODO unwrap
        let mut p = self.coeffs[ degree ];
        for i in (0..degree).rev() {
            p = p * x + self.coeffs[ i ];
        }
        p
    }

    /// Check if all the coefficients are zero
    #[inline]
    pub fn is_zero(&self) -> bool
    where
        T: Zero + PartialEq,
    {
        for i in 0..self.coeffs.len() {
            if self.coeffs[ i ] != T::zero() { return false; }
        }
        true
    }

    /// Remove leading zeros from the polynomial
    #[inline]
    pub fn trim(&mut self) 
    where
        T: Zero + PartialEq,
    {
        let mut i = self.coeffs.len() - 1;
        while self.coeffs[ i ] == T::zero() && i > 0 {
            self.coeffs.pop();
            i -= 1;
        }
    }
}

impl<T: Clone> Clone for Polynomial<T> {
    /// Clone the polynomial
    #[inline]
    fn clone(&self) -> Self {
        Self::new( self.coeffs.clone() )
    }
}

impl<T: Clone + Copy + Zero + Mul<Output = T> + Add<Output = T>> Polynomial<T> {
    /// Return the derivative of the polynomial
    #[inline]
    pub fn derivative(&self) -> Polynomial<T> {
        let mut p = Polynomial::<T>::empty();
        let degree = self.degree().unwrap(); //TODO unwrap
        p.coeffs = vec![ T::zero(); degree ];
        for i in 0..degree {
            //p.coeffs[ i ] = self.coeffs[ i + 1 ].clone() * ( i + 1 ) as f64;
            for _ in 0..=i {
                p.coeffs[ i ] = p.coeffs[ i ] + self.coeffs[ i + 1 ].clone();
            }
        }
        p
    }

    /// Return the nth derivative of the polynomial
    #[inline]
    pub fn derivative_n(&self, n: usize) -> Polynomial<T> {
        let mut p = self.clone();
        for _ in 0..n {
            p = p.derivative();
        }
        p
    }

    /// Return the nth derivative of the polynomial at a given point
    #[inline]
    pub fn derivative_at(&self, x: T, n: usize) -> T {
        let p = self.derivative_n( n );
        p.eval( x )
    }
}

impl<T: fmt::Display + Zero + std::cmp::PartialOrd + Signed> fmt::Display for Polynomial<T>
{
    /// Format the output of the polynomial
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        let degree = self.degree().unwrap(); //TODO unwrap
        for i in (0..=degree).rev() {
            let sign = if self.coeffs[ i ] < T::zero() { "-" } else { "+" };
            if i == degree && i != 1 && i != 0 {
                s.push_str( &format!( "{}x^{}", self.format_leading_coeff( i ), i ) );
                //s.push_str( &format!( "{}x^{}", self.coeffs[ i ], i ) );
            } else if i == degree && i == 1 {
                s.push_str( &format!( "{}x", self.format_leading_coeff( i ) ) );
                //s.push_str( &format!( "{}x", self.coeffs[ i ] ) );
            } else if i == degree && i == 0 {
                s.push_str( &format!( "{}", self.coeffs[ i ] ) );
            } else if i == 1 {
                s.push_str( &format!( " {} {}x", sign, self.coeffs[ i ].abs() ) );
            } else if i == 0 {
                s.push_str( &format!( " {} {}", sign, self.coeffs[ i ].abs() ) );
            } else {
                s.push_str( &format!( " {} {}x^{}", sign, self.coeffs[ i ].abs(), i ) );
            }
        }
        write!(f, "{}", s )
    }
} 

impl<T: Zero + std::cmp::PartialOrd + Signed + fmt::Display> Polynomial<T> {
    fn format_leading_coeff( &self, i: usize ) -> String {  
        let mut str = String::new();
        if self.coeffs[ i ] == T::one() {
            str.push_str( "" );
        } else if self.coeffs[ i ] == -T::one() {
            str.push_str( "-" );
        } else {
            str.push_str( &format!( "{}", self.coeffs[ i ] ) );
        }
        str
    }
}

impl<T: fmt::Debug> fmt::Debug for Polynomial<T>
{
    /// Format the debug output of the polynomial
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.coeffs )
    }
} 

impl Polynomial<f64> {
    /// Find all the roots of the polynomial with f64 coefficients
    /// refine = true will refine the roots using Laguerre's method
    #[inline]
    pub fn roots( &self, refine: bool ) -> Vector::<Cmplx> {
        let mut coeffs: Vector::<Cmplx> = Vector::new( self.coeffs.len(), Cmplx::zero() );
        for i in 0..self.coeffs.len() {
            coeffs[i] = Cmplx::new( self.coeffs[i], 0.0 ); // Convert to Complex<f64>
        }
        Polynomial::<Cmplx>::poly_solve( coeffs, refine )
    }
}

impl Polynomial<Cmplx> {
    /// Find all the roots of the polynomial with Complex<f64> coefficients
    /// refine = true will refine the roots using Laguerre's method
    #[inline]
    pub fn roots( &self, refine: bool ) -> Vector::<Cmplx> {
        let mut coeffs: Vector::<Cmplx> = Vector::new( self.coeffs.len(), Cmplx::zero() );
        for i in 0..self.coeffs.len() { coeffs[i] = self.coeffs[i]; }
        Polynomial::<Cmplx>::poly_solve( coeffs, refine )
    }

    // Solve the quadratic equation ax^2 + bx + c = 0
    fn quadratic_solve( a: Cmplx, b: Cmplx, c: Cmplx ) -> Vector::<Cmplx> {
        let mut roots = Vector::<Cmplx>::zeros( 2 );
        let discriminant: Cmplx = b * b - a * c * 4.0;
        let mut sgn: f64 = ( b.conj() * discriminant.sqrt() ).real;
        if sgn >= 0.0 { sgn = 1.0; } else { sgn = -1.0; }
        let q: Cmplx = -( b + discriminant.sqrt() * sgn ) * 0.5;
        roots[0] = q / a;
        roots[1] = c / q;
        roots
    
    }

    // Solve the cubic equation ax^3 + bx^2 + cx + d = 0
    fn cubic_solve( a: Cmplx, b: Cmplx, c: Cmplx, d: Cmplx ) -> Vector::<Cmplx> {
        let mut roots = Vector::<Cmplx>::zeros( 3 );
        let ( a2, b2, c2, d2 ) = ( a * a, b * b, c * c, d * d );
        let dis = a * b * c * d * 18. - b * b2 * d * 4. + b2 * c2 - a * c2 * c * 4. - a2 * d2 * 27.;
        let d0 = b2 - a * c * 3.;
        let d1 = b2 * b * 2. - a * b * c * 9. + a2 * d * 27.;
        if d0 == Cmplx::zero() && d1 == Cmplx::zero() { // 3 equal roots
            roots[0] = -b / ( a * 3.0 );
            roots[1] = roots[0];
            roots[2] = roots[0];
        } else {
            let sqrt = (-a * a * dis * 27.0).sqrt();
            let base = if d1 < Cmplx::zero() { d1 - sqrt } else { d1 + sqrt } / 2.0;
            let k = base.pow( &Cmplx::new( 1. / 3.0, 0.0 ) );
            roots[0] = -(b + k + d0 / k) / ( a * 3.0 );
            let u = Cmplx::new( -0.5, (3.0_f64).sqrt() / 2.0 );
            roots[1] = -(b + u * k + d0 / ( u * k ) ) / ( a * 3.0 );
            let u2 = u * u;
            roots[2] = -(b + u2 * k + d0 / ( u2 * k ) ) / ( a * 3.0 );
        }
        roots
    }

    // Find all the roots of the polynomial with Complex<f64> coefficients
    // refine = true will refine the roots using Laguerre's method
    #[inline]
    fn poly_solve( coeffs: Vector::<Cmplx>, refine: bool ) -> Vector::<Cmplx> {
        let degree = coeffs.size() - 1;
        let mut poly_roots = Vector::<Cmplx>::zeros( degree );
        if degree == 0 { panic!( "Polynomial roots error: degree must be at least one.") }
        if degree == 1 {
            poly_roots[0] = - coeffs[0] / coeffs[1];
        }
        if degree == 2 {
            let a = coeffs[2];
            let b = coeffs[1];
            let c = coeffs[0];
            poly_roots = Polynomial::quadratic_solve( a, b, c );
        }
        if degree == 3 {
            let a = coeffs[3];
            let b = coeffs[2];
            let c = coeffs[1];
            let d = coeffs[0];
            poly_roots = Polynomial::cubic_solve( a, b, c, d );
        }
        let eps = f64::EPSILON;
        let mut its: usize = 0;
        let mut a = coeffs.clone();
        if degree > 3 {
            let mut ad = coeffs.clone();
            let mut b;
            for j in (0..degree).rev() {
                let mut x = Cmplx::zero();
                let mut ad_v = Vector::<Cmplx>::zeros( j + 2 );
                for jj in 0..j+2 {
                    ad_v[jj] = ad[jj];
                }
                Self::laguer( &mut ad_v, &mut x, &mut its );
                if x.imag.abs() <= 2.0 * eps * x.real.abs() {
                    x = Cmplx::new( x.real, 0.0 );
                }
                poly_roots[j] = x;
                b = ad[ j + 1 ];
                for jj in (0..j+1).rev() {
                    let c = ad[jj];
                    ad[jj] = b;
                    b = x * b + c;
                }
            }
        }
        if refine {
            for j in 0..degree {
                Self::laguer( &mut a, &mut poly_roots[j], &mut its );
            }
        }
        poly_roots
    }

    fn laguer( a: &mut Vector::<Cmplx>, x: &mut Cmplx, iterations: &mut usize ) {
        const MR: usize = 8;
        const MT: usize = 10;
        const MAXIT: usize = MT * MR;
        const EPS: f64 = f64::EPSILON;
        let frac: [f64; MR + 1] = [ 0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0 ];
        let m = a.size() - 1;
        for iter in 1..MAXIT {
            *iterations = iter;
            let mut b = a[m];
            let mut err = b.abs();
            let mut d = Cmplx::zero();
            let mut f = Cmplx::zero();
            let abx = x.abs();
            for j in (0..m).rev() {
                f = *x * f + d;
                d = *x * d + b;
                b = *x * b + a[j];
                err = b.abs() + abx * err;
            }
            err *= EPS;
            if b.abs() <= err { return; }
            let g = d / b;
            let g2 = g * g;
            let h = g2 - ( f / b ) * 2.0;
            let sq = ( ( h * (m as f64) - g2 ) * ( m - 1 ) as f64  ).sqrt();
            let mut gp = g + sq;
            let gm = g - sq;
            let abp = gp.abs();
            let abm = gm.abs();
            if abp < abm { gp = gm; }
            let dx = if f64::max( abp, abm ) > 0.0 { 
                Cmplx::new( m as f64, 0.0 ) / gp
            } else {
                Cmplx::polar( 1.0 + abx, iter as f64 )
            };
            let x1 = *x - dx;
            if *x == x1 { return; }
            if iter % MT != 0 { *x = x1; } else { *x -= dx * frac[ iter / MT ]; }
        }
    }
}

