use core::ops::{Index, IndexMut, Add};
//use core::ops::{Neg, Add, Sub, Mul, Div};
use crate::{ Cmplx, Vector };
use crate::traits::{Zero, Number};

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
    pub fn evaluate(&self, x: T) -> T
    where
        T: Copy + core::ops::Mul<Output = T> + core::ops::Add<Output = T>,
    {
        let degree = self.degree().unwrap(); //TODO unwrap
        let mut p = self.coeffs[ degree ];
        for i in (0..degree).rev() {
            p = p * x + self.coeffs[ i ];
        }
        p
    }

    //TODO polydiv, derivative (and integral?) use Luna Polynomial class as a guide

}

impl<T: Clone> Clone for Polynomial<T> {
    /// Clone the polynomial
    #[inline]
    fn clone(&self) -> Self {
        Self::new( self.coeffs.clone() )
    }
}

impl<T: Copy + Clone + Number> Add<Polynomial<T>> for Polynomial<T> {
    type Output = Self;
    /// Add two polynomials together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        let mut sum = Polynomial::<T>::empty();
        let mut degree = match self.degree() {
            Ok( d ) => d,
            Err( _ ) => { return plus.clone(); },
        };
        let plus_degree = match plus.degree() {
            Ok( d ) => d,
            Err( _ ) => { return self.clone(); },
        };
        if degree < plus_degree {
            degree = plus_degree;
        }
        sum.coeffs = vec![ T::zero(); degree + 1 ];
        for i in 0..=degree {
            if i <= self.degree().unwrap() {
                sum.coeffs[ i ] = sum.coeffs[ i ] + self.coeffs[ i ].clone();
            }
            if i <= plus.degree().unwrap() {
                sum.coeffs[ i ] = sum.coeffs[ i ] + plus.coeffs[ i ].clone();
            }
        }
        sum
    }
}

//TODO subtraction, multiplication, division (quotient and remainder)

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

    /// Divide the polynomial by another polynomial to get a quotient and remainder
    #[inline]
    pub fn polydiv( &self, v: &Polynomial<f64> ) -> Result<( Polynomial<f64>, Polynomial<f64> ), &'static str> {
        let mut n = self.degree()?;
        let mut nv = v.degree()? as isize;
        while ( nv >= 0 ) && ( v.coeffs[ nv as usize ] == 0.0 ) { nv -= 1; }
        if nv < 0 { return Err( "Polynomial.polydiv() divide by zero polynomial" ); }
        
        let mut r = Polynomial::<f64>::new( self.coeffs.clone() );
        let mut q = Polynomial::<f64>::empty();
        q.coeffs = vec![ 0.0; n - nv as usize + 1 ];

        for k in (0..=n-nv as usize).rev() {
            q.coeffs[ k ] = r.coeffs[ nv as usize + k ] / v.coeffs[ nv as usize ];
            for j in (k..=nv as usize).rev() {
                r.coeffs[ j ] -= q.coeffs[ k ] * v.coeffs[ j - k ];
            }
        }
        for j in (nv as usize)..=n {
            r.coeffs[ j ] = 0.0;
        }
        // Remove leading zeros from r and q polynomials
        let mut r_val = r.coeffs[ n ];
        while r_val.abs() == 0.0 && n > 0 {
            r.coeffs.pop();
            n -= 1;
            r_val = r.coeffs[ n ];
        }
        n = self.degree()? - nv as usize;
        let mut q_val = q.coeffs[ n ];
        while q_val.abs() == 0.0 && n > 0 {
            q.coeffs.pop();
            n -= 1;
            q_val = q.coeffs[ n ];
        }
        Ok( ( q, r ) )
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

impl<T> Index<usize> for Polynomial<T> {
    type Output = T;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: usize ) -> &'a T {
        if index >= self.coeffs.len() { panic!( "Index out of bounds" ); }
        &self.coeffs[ index ]
    }
}

impl<T> IndexMut<usize> for Polynomial<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: usize ) -> &mut T {
        if index >= self.coeffs.len() { panic!( "Index out of bounds" ); }
        &mut self.coeffs[ index ] 
    }
}