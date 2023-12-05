use std::fmt;
use core::ops::{ 
    Index, IndexMut, Neg, Add, Sub, Mul, Div, 
    AddAssign, SubAssign, MulAssign, DivAssign
};
use crate::Matrix;
use crate::Complex;
use crate::vector::Vector;
use crate::traits::{Zero, Number, Signed};

pub struct Tridiagonal<T> {
    sub: Vector<T>,       // Subdiagonal of the matrix (size n-1)
    main: Vector<T>,       // Main diagonal of the matrix (size n)
    sup: Vector<T>,       // Superdiagonal of the matrix (size n-1)
    n: usize            // Size of the matrix (n x n) and main diagonal
}

impl<T> Tridiagonal<T> {
    /// Create a new tridiagonal matrix of unspecified size
    #[inline]
    pub const fn empty() -> Self {
        let sub = Vector::<T>::empty();
        let main = Vector::<T>::empty();
        let sup = Vector::<T>::empty();
        let n = 0;
        Tridiagonal { sub, main, sup, n }
    }

    /// Create a new tridiagonal matrix from three vectors
    /// The vectors must be of size n-1, n and n-1 respectively.
    #[inline]
    pub fn with_vectors( sub: Vector<T>, main: Vector<T>, sup: Vector<T> ) -> Self {
        let n = main.size();
        if sub.size() != n - 1 || sup.size() != n - 1 { 
            panic!("Tridiagonal error: invalid vector sizes."); 
        }
        Tridiagonal { sub, main, sup, n }
    }

    /// Create a new tridiagonal matrix from three vecs 
    /// The vecs must be of size n-1, n and n-1 respectively.
    #[inline]
    pub fn with_vecs( sub: Vec<T>, main: Vec<T>, sup: Vec<T> ) -> Self {
        let n = main.len();
        if sub.len() != n - 1 || sup.len() != n - 1 { 
            panic!("Tridiagonal error: invalid vector sizes."); 
        }
        Tridiagonal { 
            sub: Vector::create( sub ), 
            main: Vector::create( main ), 
            sup: Vector::create( sup ), n 
        }
    }

    /// Return the size of the matrix
    #[inline]
    pub const fn size( &self ) -> usize {
        self.n
    }

    /// Return the subdiagonal of the matrix
    #[inline]
    pub fn subdiagonal( &self ) -> &Vector<T> {
        &self.sub
    }

    /// Return the main diagonal of the matrix
    #[inline]
    pub fn maindiagonal( &self ) -> &Vector<T> {
        &self.main
    }

    /// Return the superdiagonal of the matrix
    #[inline]
    pub fn superdiagonal( &self ) -> &Vector<T> {
        &self.sup
    }
}

impl<T: Clone + Copy + Zero + Number> Tridiagonal<T> {
    /// Create a new tridiagonal matrix of specified size
    #[inline]
    pub fn new( n: usize ) -> Self {
        let sub = Vector::<T>::new( n - 1, T::zero());
        let main = Vector::<T>::new( n, T::zero());
        let sup = Vector::<T>::new( n - 1, T::zero());
        Tridiagonal { sub, main, sup, n }
    }

    /// Create a new tridiagonal matrix from three elements
    #[inline]
    pub fn with_elements( sub: T, main: T, sup: T, n: usize ) -> Self {
        let sub = Vector::<T>::new( n - 1, sub);
        let main = Vector::<T>::new( n, main);
        let sup = Vector::<T>::new( n - 1, sup);
        Tridiagonal { sub, main, sup, n }
    }

    /// Resize the matrix
    #[inline]
    pub fn resize( &mut self, n: usize ) {
        //TODO maybe keep the old values?
        self.sub = Vector::<T>::new(n - 1, T::zero());
        self.main = Vector::<T>::new(n, T::zero());
        self.sup = Vector::<T>::new(n - 1, T::zero());
        self.n = n;
    }

    /// Transpose the matrix in place
    #[inline]
    pub fn transpose_in_place( &mut self ) {
        let temp = self.sub.clone();
        self.sub = self.sup.clone();
        self.sup = temp;
    }

    /// Transpose the matrix
    #[inline]
    pub fn transpose( &self ) -> Self {
        let mut temp = self.clone();
        temp.transpose_in_place();
        temp
    }

    /// Retrun the determinant of the matrix
    #[inline]
    pub fn det( &self ) -> T {
        let mut f = Vector::new( self.n + 1, T::zero() );
        f[ 0 ] = T::one();
        f[ 1 ] = self.main[ 0 ] * f[ 0 ];
        for j in 2..self.n + 1 {
            f[ j ] = self.main[ j - 1 ] * f[ j - 1 ] 
                   - self.sub[ j - 2 ] * self.sup[ j - 2 ] * f[ j - 2 ];
        }
        f[ self.n ]
    }

    /// Convert the tridiagonal matrix to a dense matrix
    #[inline]
    pub fn convert( &self ) -> Matrix<T> {
        let mut dense = Matrix::<T>::new( self.n, self.n, T::zero() );
        if self.n == 0 { panic!("Tridiagonal error: zero size matrix."); }
        if self.n == 1 {
            dense[0][0] = self.main[0];
        } else {
            dense[0][0] = self.main[0];
            dense[0][1] = self.sup[0];
            for i in 1..self.n - 1 {
                dense[i][i - 1] = self.sub[i - 1];
                dense[i][i] = self.main[i];
                dense[i][i + 1] = self.sup[i];
            }
            dense[self.n - 1][self.n - 2] = self.sub[self.n - 2];
            dense[self.n - 1][self.n - 1] = self.main[self.n - 1];
        }
        dense
    }

    /// Solve the system of equations Tx=r where r is a specified Vector
    #[inline]
    pub fn solve( &self, r: &Vector<T> ) -> Vector<T> {
        if self.n != r.size() { panic!( "Tridiagonal error: matrix and vector sizes do not agree." ); }
        let mut u = Vector::<T>::new( self.n, T::zero() );
        let mut a_temp = self.sub.clone();
        let mut c_temp = self.sup.clone();
        a_temp.push_front( T::zero() );
        c_temp.push( T::zero() );
        let mut beta = self.main[0];
        let mut gamma = Vector::<T>::new( self.n, T::zero() );
        if self.main[0] == T::zero() { panic!( "Tridiagonal error: zero on leading diagonal." ); }
        u[0] = r[0] / beta;
        for j in 1..self.n {
            gamma[j] = c_temp[j - 1] / beta;
            beta = self.main[j] - a_temp[j] * gamma[j];
            if beta == T::zero() { panic!( "Tridiagonal error: zero pivot." ); }
            u[j] = ( r[j] - a_temp[j] * u[j - 1] ) / beta;
        }
        for j in (0..self.n - 1).rev() {
            let temp = gamma[j + 1] * u[j + 1];
            u[j] -= temp;
        }
        u
    }
}

impl<T: Clone + Signed> Tridiagonal::<Complex::<T>> {
    /// Return the conjugate of the tridiagonal matrix
    #[inline]
    pub fn conj( &self ) -> Self {
        let sub = self.sub.conj();
        let main = self.main.conj();
        let sup = self.sup.conj();
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: fmt::Debug> fmt::Debug for Tridiagonal<T>
{
    /// Format the debug output of the tridiagonal matrix
    fn fmt( &self, f: &mut fmt::Formatter<'_> ) -> fmt::Result {
        write!(f, "n: {}\n", self.n)?;
        write!(f, "sub: {:?}\n", self.sub)?;
        write!(f, "main: {:?}\n", self.main)?;
        write!(f, "super: {:?}\n", self.sup)?;
        write!(f, "\n")
    }
} 

impl<T: fmt::Display> fmt::Display for Tridiagonal<T> 
{
    /// Format the output of the tridiagonal matrix
    fn fmt( &self, f: &mut fmt::Formatter<'_> ) -> fmt::Result {
        if self.n == 0 { return write!(f, "Empty tridiagonal matrix\n"); }
        for j in 0..self.n {
            let mut tabs = String::new();
            let mut tabs_right = String::new();
            for _i in 1..j { tabs += "*\t"; }
            for _i in 1..self.n - j - 1 { tabs_right += "*\t"; }
            if j == 0 {
                let string = format!( "{}\t{}\t{}", self.main[0], self.sup[0], tabs_right );
                writeln!(f,"{}", string )?;
            } else if j == self.n - 1 {
                let string = format!( "{}{}\t{}", tabs, self.sub[self.n - 2], self.main[self.n - 1] );
                writeln!(f,"{}", string )?;
            } else {
                let string = format!( "{}{}\t{}\t{}\t{}", 
                    tabs, self.sub[j - 1], self.main[j], self.sup[j], tabs_right );
                writeln!(f,"{}", string )?;
            }
        }
        write!(f, "\n")
    }
}

impl<T> Index<(usize, usize)> for Tridiagonal<T> {
    type Output = T;
    /// Indexing operator [(i,j)] (read only)
    #[inline]
    fn index<'a>(&'a self, index: ( usize, usize ) ) -> &'a T {
        let ( i, j ) = index;
        if i >= self.n || j >= self.n { panic!("Tridiagonal error: index out of bounds."); }
        if i == j { return &self.main[i]; }
        if i == j + 1 { return &self.sub[j]; }
        if i + 1 == j { return &self.sup[i]; }
        panic!("Tridiagonal error: index out of bounds.");
    }
}

impl<T> IndexMut<(usize, usize)> for Tridiagonal<T> {
    /// Indexing operator [(i,j)] (read/write)
    #[inline]
    fn index_mut(&mut self, index: ( usize, usize ) ) -> &mut T {
        let ( i, j ) = index;
        if i >= self.n || j >= self.n { panic!("Tridiagonal error: index out of bounds."); }
        if i == j { return &mut self.main[i]; }
        if i == j + 1 { return &mut self.sub[j]; }
        if i + 1 == j { return &mut self.sup[i]; }
        panic!("Tridiagonal error: index out of bounds."); 
    }
}

impl<T: Clone> Clone for Tridiagonal<T> {
    /// Clone the tridiagonal matrix
    #[inline]
    fn clone(&self) -> Self {
        let sub = self.sub.clone();
        let main = self.main.clone();
        let sup = self.sup.clone();
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Neg<Output = T>> Neg for Tridiagonal<T> {
    type Output = Self;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        let sub = -self.sub;
        let main = -self.main;
        let sup = -self.sup;
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Number> Add<Tridiagonal<T>> for Tridiagonal<T> {
    type Output = Self;
    /// Add the elements of two tridiagonal matrices together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        if self.size() != plus.size() { panic!( "Tridiagonal matrix sizes do not agree (+)." ); }
        let sub = self.sub + plus.sub;
        let main = self.main + plus.main;
        let sup = self.sup + plus.sup;
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Number> Sub<Tridiagonal<T>> for Tridiagonal<T> {
    type Output = Self;
    /// Subtract the elements of one tridiagonal matrix from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        if self.size() != minus.size() { panic!( "Tridiagonal matrix sizes do not agree (-)." ); }
        let sub = self.sub - minus.sub;
        let main = self.main - minus.main;
        let sup = self.sup - minus.sup;
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Number> Mul<T> for Tridiagonal<T> {
    type Output = Self;
    /// Multiply a tridiagonal matrix by a scalar (tridiagonal * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        let sub = self.sub * scalar.clone();
        let main = self.main * scalar.clone();
        let sup = self.sup * scalar.clone();
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl Mul<Tridiagonal<f64>> for f64 {
    type Output = Tridiagonal<f64>;
    /// Allow multiplication on the left by f64 (f64 * tridiagonal matrix)
    #[inline]
    fn mul(self, tridiagonal: Tridiagonal<f64>) -> Self::Output {
        let sub = self * tridiagonal.sub;
        let main = self * tridiagonal.main;
        let sup = self * tridiagonal.sup;
        let n = tridiagonal.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Number> Div<T> for Tridiagonal<T> {
    type Output = Self;
    /// Divide a vector by a scalar (vector / scalar)
    fn div(self, scalar: T) -> Self::Output {
        let sub = self.sub / scalar.clone();
        let main = self.main / scalar.clone();
        let sup = self.sup / scalar.clone();
        let n = self.n;
        Tridiagonal { sub, main, sup, n }
    }
}

impl<T: Clone + Number> AddAssign<T> for Tridiagonal<T> {
    /// Add the same value to every element in a mutable tridiagonal matrix
    #[inline]
    fn add_assign(&mut self, plus: T) {
        self.sub += plus.clone();
        self.main += plus.clone();
        self.sup += plus.clone();
    }
}

impl<T: Clone + Number> SubAssign<T> for Tridiagonal<T> {
    /// Subtract the same value from every element in a mutable tridiagonal matrix 
    fn sub_assign(&mut self, rhs: T) {
        self.sub -= rhs.clone();
        self.main -= rhs.clone();
        self.sup -= rhs.clone();
    }
} 

impl<T: Clone + Number> MulAssign<T> for Tridiagonal<T> {
    /// Multiply every element in a mutable tridiagonal matrix by a scalar 
    fn mul_assign(&mut self, rhs: T) {
        self.sub *= rhs.clone();
        self.main *= rhs.clone();
        self.sup *= rhs.clone();
    }
} 

impl<T: Clone + Number> DivAssign<T> for Tridiagonal<T> {
    /// Divide every element in a mutable tridiagonal matrix by a scalar 
    fn div_assign(&mut self, rhs: T) {
        self.sub /= rhs.clone();
        self.main /= rhs.clone();
        self.sup /= rhs.clone();
    }
}

impl<T: Clone + Copy + Number> Mul<Vector<T>> for Tridiagonal<T> {
    type Output = Vector<T>;
    /// Multiply a tridiagonal matrix with a (column) vector ( tridiagonal matrix * vector )
    #[inline]
    fn mul(self, vec: Vector<T> ) -> Vector<T> {
        if self.size() != vec.size() { 
            panic!( "Tridiagonal matrix and vector sizes do not agree (*)." ); 
        }
        let mut result = Vector::<T>::new( self.size(), T::zero() );
        result[ 0 ] = self.main[ 0 ] * vec[ 0 ] + self.sup[ 0 ] * vec[ 1 ];
        for i in 1..self.size() - 1 {
            result[ i ] = self.sub[ i - 1 ] * vec[ i - 1 ] + self.main[ i ] * vec[ i ]
                        + self.sup[ i ] * vec[ i + 1 ];
        }
        result[ self.n - 1 ] = self.sub[ self.n - 2 ] * vec[ self.n - 2 ]  
                             + self.main[ self.n - 1 ] * vec[ self.n - 1 ];
        result
    }
}

// Non-consuming vector multiplication operator
impl<T: Clone + Copy + Number> Mul<&Vector<T>> for &Tridiagonal<T> {
    type Output = Vector<T>;
    /// Multiply a tridiagonal matrix with a (column) vector ( tridiagonal matrix * vector )
    #[inline]
    fn mul(self, vec: &Vector<T>) -> Vector<T> {
        if self.size() != vec.size() { 
            panic!( "Tridiagonal matrix and vector sizes do not agree (*)." ); 
        }
        let mut result = Vector::<T>::new( self.size(), T::zero() );
        result[ 0 ] = self.main[ 0 ] * vec[ 0 ] + self.sup[ 0 ] * vec[ 1 ];
        for i in 1..self.size() - 1 {
            result[ i ] = self.sub[ i - 1 ] * vec[ i - 1 ] + self.main[ i ] * vec[ i ]
                        + self.sup[ i ] * vec[ i + 1 ];
        }
        result[ self.n - 1 ] = self.sub[ self.n - 2 ] * vec[ self.n - 2 ]  
                             + self.main[ self.n - 1 ] * vec[ self.n - 1 ];
        result
        
    }
}

