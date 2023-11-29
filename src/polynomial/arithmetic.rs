use crate::polynomial::Polynomial;
use core::ops::{Index, IndexMut, Add, Sub, Neg, Mul};
use crate::traits::{Number, Signed};

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
        if degree < plus_degree { degree = plus_degree; }
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

impl<T: Clone + Signed> Neg for Polynomial<T> {
    type Output = Self;
    /// Return the unary negation ( unary - )
    #[inline]
    fn neg(self) -> Self::Output {
        let mut neg = Polynomial::<T>::empty();
        neg.coeffs = self.coeffs.iter().map( |x| -x.clone() ).collect();
        neg
    }
}

impl<T: Copy + Clone + Number + Signed> Sub<Polynomial<T>> for Polynomial<T> {
    type Output = Self;
    /// Subtract one polynomial from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        let mut diff = Polynomial::<T>::empty();
        let mut degree = match self.degree() {
            Ok( d ) => d,
            Err( _ ) => { return - minus.clone(); },
        };
        let minus_degree = match minus.degree() {
            Ok( d ) => d,
            Err( _ ) => { return self.clone(); },
        };
        if degree < minus_degree { degree = minus_degree; }
        diff.coeffs = vec![ T::zero(); degree + 1 ];
        for i in 0..=degree {
            if i <= self.degree().unwrap() {
                diff.coeffs[ i ] = diff.coeffs[ i ] + self.coeffs[ i ].clone();
            }
            if i <= minus.degree().unwrap() {
                diff.coeffs[ i ] = diff.coeffs[ i ] - minus.coeffs[ i ].clone();
            }
        }
        diff
    }
}

impl<T: Copy + Clone + Number> Mul<Polynomial<T>> for Polynomial<T> {
    type Output = Self;
    /// Multiply two polynomials together ( binary * )
    #[inline]
    fn mul(self, times: Self) -> Self::Output {
        let mut product = Polynomial::<T>::empty();
        let mut degree = match self.degree() {
            Ok( d ) => d,
            Err( _ ) => { return Polynomial::<T>::empty(); },
        };
        let times_degree = match times.degree() {
            Ok( d ) => d,
            Err( _ ) => { return Polynomial::<T>::empty(); },
        };
        degree += times_degree;
        product.coeffs = vec![ T::zero(); degree + 1 ];
        for i in 0..=self.degree().unwrap() {
            for j in 0..=times.degree().unwrap() {
                product.coeffs[ i + j ] = product.coeffs[ i + j ] + self.coeffs[ i ].clone() * times.coeffs[ j ].clone();
            }
        }
        product
    }
}

impl<T: Copy + Clone + Number> Mul<T> for Polynomial<T> {
    type Output = Self;
    /// Multiply a polynomial by a scalar ( binary * )
    #[inline]
    fn mul(self, times: T) -> Self::Output {
        let mut product = Polynomial::<T>::empty();
        product.coeffs = self.coeffs.iter().map( |x| x.clone() * times.clone() ).collect();
        product
    }
}

impl<T: Copy + Clone + Number + Signed + std::fmt::Debug> Polynomial<T> {
    /// Divide the polynomial by another polynomial to get a quotient and remainder
    #[inline]
    pub fn polydiv( &self, v: &Polynomial<T> ) -> Result<( Polynomial<T>, Polynomial<T> ), &'static str> {
        if v.coeffs.len() == 0 { return Err( "Polynomial.polydiv() divide by zero polynomial" ); }
        if v.is_zero() { return Err( "Polynomial.polydiv() divide by zero polynomial" ); }
        let mut q = Polynomial::<T>::empty();
        let mut r = self.clone();
        const MAX: usize = 1000; // TODO make this a parameter in the struct?
        let mut count = 0;
        while !r.is_zero() && r.degree()? >= v.degree()? {
            let mut t = Polynomial::<T>::empty();
            t.coeffs = vec![ T::zero(); r.degree()? - v.degree()? + 1 ];
            t.coeffs[ r.degree()? - v.degree()? ] = r.coeffs[ r.degree()? ] / v.coeffs[ v.degree()? ];
            q = q + t.clone();
            r = r - ( t * v.clone() );
            r.trim();
            q.trim();
            count += 1;
            if count > MAX { return Err( "Polynomial.polydiv() exceeded maximum iterations" ); }
        }
        Ok( ( q ,r ) )
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