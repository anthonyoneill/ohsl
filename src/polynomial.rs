use core::ops::{Index, IndexMut};
//use core::ops::{Neg, Add, Sub, Mul, Div};
use crate::{ Cmplx, Vector };
use crate::traits::Zero;

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

    /// Return the degree of the polynomial
    #[inline]
    pub fn degree(&self) -> usize {
        self.coeffs.len() - 1
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
        let mut p = self.coeffs[ self.degree() ];
        for i in (0..self.degree()).rev() {
            p = p * x + self.coeffs[ i ];
        }
        p
    }

    

    /// Divide the polynomial by another polynomial to get a quotient and remainder
    #[inline]
    pub fn polydiv( &self, _v: &Polynomial<T> ) -> ( Polynomial<T>, Polynomial<T> ) {
        //TODO
        ( Polynomial::empty(), Polynomial::empty() )
    }

    //TODO polydiv, solve/roots, derivative (and integral?) use Luna Polynomial class as a guide

}

impl Polynomial<f64> {
    /// Find all the roots of the polynomial with f64 coefficients
    /// refine = true will refine the roots using Laguerre's method
    #[inline]
    pub fn roots( &self, _refine: bool ) -> Vector::<Cmplx> {
        let mut poly_roots = Vector::<Cmplx>::zeros( self.degree() );
        if self.degree() == 0 { panic!( "Polynomial roots error: degree must be at least one.") }
        if self.degree() == 1 {
            poly_roots[0] = Cmplx::new( - self.coeffs[0] / self.coeffs[1], 0.0 );
            //return poly_roots;
        }
        if self.degree() == 2 {
            let a = Cmplx::new( self.coeffs[2], 0.0 );
            let b = Cmplx::new( self.coeffs[1], 0.0 );
            let c = Cmplx::new( self.coeffs[0], 0.0 );
            poly_roots = Polynomial::quadratic_solve( a, b, c );
            //return poly_roots;
        }
        if self.degree() == 3 {
            let a = Cmplx::new( self.coeffs[3], 0.0 );
            let b = Cmplx::new( self.coeffs[2], 0.0 );
            let c = Cmplx::new( self.coeffs[1], 0.0 );
            let d = Cmplx::new( self.coeffs[0], 0.0 );
            println!( "a = {}, b = {}, c = {}, d = {}", a, b, c, d );
            poly_roots = Polynomial::cubic_solve( a, b, c, d );
            //return poly_roots;
        }


        poly_roots
    }

}

impl Polynomial<Cmplx> {
    // Solve the quadratic equation ax^2 + bx + c = 0
    fn quadratic_solve( a: Cmplx, b: Cmplx, c: Cmplx ) -> Vector<Cmplx> {
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
    fn cubic_solve( a: Cmplx, b: Cmplx, c: Cmplx, d: Cmplx ) -> Vector<Cmplx> {
        let mut roots = Vector::<Cmplx>::zeros( 3 );
        let disc =  a * b * c * d * 18.0 - b * b * b * d * 4.0 + b * b * c * c - a * c * c * c * 4.0 - a * a * d * d * 27.0;
        let d0 = b * b - a * c * 3.0;
        let d1 = b * b * b * 2.0 - a * b * c * 9.0 + a * a * d * 27.0;
        if d0 == Cmplx::zero() && d1 == Cmplx::zero() { // 3 equal roots
            roots[0] = -b / ( a * 3.0 );
            roots[1] = roots[0];
            roots[2] = roots[0];
        } else {
            let sqrt = (-a * a * disc * 27.0).sqrt();
            let k = ( if d1 < Cmplx::zero() { d1 - sqrt } else { d1 + sqrt } / 2.0).pow( &Cmplx::new( 1. / 3.0, 0.0 ) );
            roots[0] = -(b + k + d0 / k) / ( a * 3.0 );
            let u = Cmplx::new( -0.5, (3.0_f64).sqrt() / 2.0 );
            roots[1] = -(b + u * k + d0 / ( u * k ) ) / ( a * 3.0 );
            let u2 = u * u;
            roots[2] = -(b + u2 * k + d0 / ( u2 * k ) ) / ( a * 3.0 );
        }
        roots
    }
}

//TODO impl Polynomial<Cmplx> 

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