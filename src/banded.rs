use std::fmt;
use std::ops::{Index, IndexMut, Neg, Add, Sub, Mul, Div};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use crate::{Matrix, Number};
use crate::traits::Signed;

//TODO n = size, m1 = below, m2 = above

#[derive(PartialEq, Clone)]
pub struct Banded<T> {
    n: usize,           // Size of the matrix (n x n) and main diagonal
    m1: usize,          // Number of bands below the main diagonal
    m2: usize,          // Number of bands above the main diagonal
    compact: Matrix<T>, // Compact representation of the banded matrix
}

impl<T> Banded<T> {
    /// Create a new banded matrix of unspecified size
    #[inline]
    pub fn empty() -> Self {
        Banded {
            n: 0,
            m1: 0,
            m2: 0,
            compact: Matrix::empty(),
        }
    }

    /// Return the size of the matrix
    #[inline]
    pub fn size(&self) -> usize {
        self.n
    }

    /// Return the number of bands below the main diagonal
    #[inline]
    pub fn size_below(&self) -> usize {
        self.m1
    }

    /// Return the number of bands above the main diagonal
    #[inline]
    pub fn size_above(&self) -> usize {
        self.m2
    }

    /// Return the compact representation of the banded matrix
    #[inline]
    pub fn compact(&self) -> &Matrix<T> {
        &self.compact
    }
}

impl<T: Clone + Copy + Number> Banded<T> {
    /// Create a new banded matrix of specified size and fill it with a constant value
    #[inline]
    pub fn new( n: usize, m1: usize, m2: usize, value: T ) -> Self {
        Banded {
            n,
            m1,
            m2,
            compact: Matrix::new( n, m1 + m2 + 1, value ),
        }
    }

    /// Fill the entire matrix with a constant value
    #[inline]
    pub fn fill(&mut self, value: T ) {
        self.compact.fill( value );
    }

    /// Resize the matrix to a new size
    #[inline]
    pub fn resize(&mut self, n: usize, m1: usize, m2: usize) {
        self.n = n;
        self.m1 = m1;
        self.m2 = m2;
        self.compact.resize( n, m1 + m2 + 1 );
    }

    /// Fill a particular band of the matrix with a constant value
    /// 0 = main diagonal, positive = above the main diagonal, negative = below the main diagonal
    #[inline]
    pub fn fill_band(&mut self, band: isize, value: T ) {
        if band < - (self.m1 as isize) || band > self.m2 as isize { 
            panic!("Banded error: band not in matrix."); 
        }
        self.compact.fill_col( (self.m1 as isize + band) as usize, value );
    }
}

impl<T> Index<(usize, usize)> for Banded<T> {
    type Output = T;
    /// Indexing operator [(i,j)] (read only)
    #[inline]
    fn index<'a>(&'a self, index: ( usize, usize ) ) -> &'a T {
        let ( i, j ) = index;
        if j > i + self.m2 || i > j + self.m1 { panic!("Banded error: index not in a band."); }
        &self.compact[ i ][ self.m1 + j - i ]
    }
}

impl<T> IndexMut<(usize, usize)> for Banded<T> {
    /// Indexing operator [(i,j)] (read/write)
    #[inline]
    fn index_mut(&mut self, index: ( usize, usize ) ) -> &mut T {
        let ( i, j ) = index;
        if j > i + self.m2 || i > j + self.m1 { panic!("Banded error: index not in a band."); }
        &mut self.compact[ i ][ self.m1 + j - i ]
    }
}

impl<T: fmt::Debug> fmt::Debug for Banded<T> {
    /// Format the debug output of the tridiagonal matrix
    fn fmt( &self, f: &mut fmt::Formatter<'_> ) -> fmt::Result {
        write!(f, "n: {}\n", self.n)?;
        write!(f, "m1: {:?}\n", self.m1)?;
        write!(f, "m2: {:?}\n", self.m2)?;
        write!(f, "compact:\n{:?}\n", self.compact)
    }
} 

impl<T: fmt::Display> fmt::Display for Banded<T> {
    /// Format the display output of the tridiagonal matrix
    fn fmt( &self, f: &mut fmt::Formatter<'_> ) -> fmt::Result {
        if self.n == 0 { return write!(f, "Empty tridiagonal matrix\n"); }
        for i in 0..self.n {
            for j in 0..self.n {
                if j > i + self.m2 || i > j + self.m1 { write!(f, "\t{}", "*")?; }
                else { write!(f, "\t{}", self[(i,j)])?; }
            }
            write!(f, "\n")?;
        }
        write!(f, "\n")
    }
}

/* --- Operators --- */

// Non-consuming negation
impl<T: Copy + Signed> Neg for &Banded<T> {
    type Output = Banded<T>;
    /// Return the unary negation ( unary - )
    #[inline]
    fn neg(self) -> Self::Output {
        Banded {
            n: self.n,
            m1: self.m1,
            m2: self.m2,
            compact: -&self.compact,
        }
    }
}

// Consuming negation
impl<T: Copy + Signed> Neg for Banded<T> {
    type Output = Self;
    /// Return the unary negation ( unary - )
    #[inline]
    fn neg(self) -> Self::Output {
        -&self
    }
}

// Non-consuming addition
impl<T: Copy + Clone + Number> Add<&Banded<T>> for &Banded<T> {
    type Output = Banded<T>;
    /// Add two banded matrices together ( binary + )
    #[inline]
    fn add(self, plus: &Banded<T>) -> Self::Output {
        if self.n != plus.n { panic!( "Banded matrix size dimensions do not agree (+)." ); }
        if self.m1 != plus.m1 { panic!( "Banded matrix m1 dimensions do not agree (+)." ); }
        if self.m2 != plus.m2 { panic!( "Banded matrix m2 dimensions do not agree (+)." ); }
        Banded {
            n: self.n,
            m1: self.m1,
            m2: self.m2,
            compact: &self.compact + &plus.compact,
        }
    }
}

// Consuming addition
impl<T: Copy + Clone + Number> Add<Banded<T>> for Banded<T> {
    type Output = Self;
    /// Add two banded matrices together ( binary + )
    #[inline]
    fn add(self, plus: Banded<T>) -> Self::Output {
        &self + &plus
    }
}

// Non-consuming subtraction
impl<T: Copy + Clone + Number> Sub<&Banded<T>> for &Banded<T> {
    type Output = Banded<T>;
    /// Subtract two banded matrices together ( binary - )
    #[inline]
    fn sub(self, minus: &Banded<T>) -> Self::Output {
        if self.n != minus.n { panic!( "Banded matrix size dimensions do not agree (-)." ); }
        if self.m1 != minus.m1 { panic!( "Banded matrix m1 dimensions do not agree (-)." ); }
        if self.m2 != minus.m2 { panic!( "Banded matrix m2 dimensions do not agree (-)." ); }
        Banded {
            n: self.n,
            m1: self.m1,
            m2: self.m2,
            compact: &self.compact - &minus.compact,
        }
    }
}

// Consuming subtraction
impl<T: Copy + Clone + Number> Sub<Banded<T>> for Banded<T> {
    type Output = Self;
    /// Subtract two banded matrices together ( binary - )
    #[inline]
    fn sub(self, minus: Banded<T>) -> Self::Output {
        &self - &minus
    }
}

// Non-consuming scalar multiplication
impl<T: Copy + Clone + Number> Mul<T> for &Banded<T> {
    type Output = Banded<T>;
    /// Multiply a banded matrix by a scalar ( binary * )
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        Banded {
            n: self.n,
            m1: self.m1,
            m2: self.m2,
            compact: &self.compact * scalar,
        }
    }
}

// Consuming scalar multiplication
impl<T: Copy + Clone + Number> Mul<T> for Banded<T> {
    type Output = Self;
    /// Multiply a banded matrix by a scalar ( binary * )
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        &self * scalar
    }
}

// Non-consuming scalar division
impl<T: Copy + Clone + Number> Div<T> for &Banded<T> {
    type Output = Banded<T>;
    /// Divide a banded matrix by a scalar ( binary / )
    #[inline]
    fn div(self, scalar: T) -> Self::Output {
        Banded {
            n: self.n,
            m1: self.m1,
            m2: self.m2,
            compact: &self.compact / scalar,
        }
    }
}

// Consuming scalar division
impl<T: Copy + Clone + Number> Div<T> for Banded<T> {
    type Output = Self;
    /// Divide a banded matrix by a scalar ( binary / )
    #[inline]
    fn div(self, scalar: T) -> Self::Output {
        &self / scalar
    }
}

// Non-consuming addition assignment
impl<T: Copy + Clone + Number> AddAssign<&Banded<T>> for Banded<T> {
    /// Add a banded matrix to another ( += )
    #[inline]
    fn add_assign(&mut self, plus: &Banded<T>) {
        if self.n != plus.n { panic!( "Banded matrix size dimensions do not agree (+=)." ); }
        if self.m1 != plus.m1 { panic!( "Banded matrix m1 dimensions do not agree (+=)." ); }
        if self.m2 != plus.m2 { panic!( "Banded matrix m2 dimensions do not agree (+=)." ); }
        self.compact += &plus.compact;
    }
}

// Consuming addition assignment
impl<T: Copy + Clone + Number> AddAssign<Banded<T>> for Banded<T> {
    /// Add a banded matrix to another ( += )
    #[inline]
    fn add_assign(&mut self, plus: Banded<T>) {
        *self += &plus;
    }
}

// Non-consuming subtraction assignment
impl<T: Copy + Clone + Number> SubAssign<&Banded<T>> for Banded<T> {
    /// Subtract a banded matrix from another ( -= )
    #[inline]
    fn sub_assign(&mut self, minus: &Banded<T>) {
        if self.n != minus.n { panic!( "Banded matrix size dimensions do not agree (-=)." ); }
        if self.m1 != minus.m1 { panic!( "Banded matrix m1 dimensions do not agree (-=)." ); }
        if self.m2 != minus.m2 { panic!( "Banded matrix m2 dimensions do not agree (-=)." ); }
        self.compact -= &minus.compact;
    }
}

// Consuming subtraction assignment
impl<T: Copy + Clone + Number> SubAssign<Banded<T>> for Banded<T> {
    /// Subtract a banded matrix from another ( -= )
    #[inline]
    fn sub_assign(&mut self, minus: Banded<T>) {
        *self -= &minus;
    }
}

// Scalar multiplication assignment
impl<T: Copy + Clone + Number> MulAssign<T> for Banded<T> {
    /// Multiply a banded matrix by a scalar and assign the result ( *= )
    #[inline]
    fn mul_assign(&mut self, scalar: T) {
        self.compact *= scalar;
    }
}

// Scalar division assignment
impl<T: Copy + Clone + Number> DivAssign<T> for Banded<T> {
    /// Divide a banded matrix by a scalar and assign the result ( /= )
    #[inline]
    fn div_assign(&mut self, scalar: T) {
        self.compact /= scalar;
    }
}

// Constant addition assignment
impl<T: Copy + Clone + Number> AddAssign<T> for Banded<T> {
    /// Add a constant to a banded matrix and assign the result ( += )
    #[inline]
    fn add_assign(&mut self, constant: T) {
        self.compact += constant;
    }
}

// Constant subtraction assignment
impl<T: Copy + Clone + Number> SubAssign<T> for Banded<T> {
    /// Subtract a constant from a banded matrix and assign the result ( -= )
    #[inline]
    fn sub_assign(&mut self, constant: T) {
        self.compact -= constant;
    }
}

//TODO operators (banded matrix * vector ), solve, determinant