use core::ops::{Index, IndexMut, Neg, Add, Sub, Mul, Div};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::{fmt, fs::File, io::Write, cmp::Ordering};
use rand::Rng;

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Complex;

pub struct Vector<T> {
    pub vec: Vec<T>,
    size: usize,
}

pub type Vec64 = Vector<f64>;

impl<T> Vector<T> {
    /// Create a new vector of unspecified size
    #[inline]
    pub const fn empty() -> Self {
        let vec = Vec::<T>::new();
        let size = 0;
        Vector { vec, size }
    }

    /// Create a vector from an std::vec::Vec
    #[inline]
    pub fn create( vec: Vec<T> ) -> Self {
        let size = vec.len();
        Vector { vec, size }
    }

    /// Return the size of the vector 
    #[inline]
    pub fn size(&self) -> usize {
        self.size
    }

    /// Remove all of the element from the vector
    #[inline]
    pub fn clear(&mut self) {
        self.vec.clear();
        self.size = 0;
    }

    /// Swap elements two elements in the vector 
    #[inline]
    pub fn swap(&mut self, i: usize, j: usize) {
        self.vec.swap( i, j );
    }

    /// Push a new element onto the end of the vector 
    #[inline]
    pub fn push(&mut self, elem: T ) {
        self.vec.push( elem );
        self.size += 1;
    }

    /// Insert a new element into the Vector at a specified position 
    #[inline]
    pub fn insert(&mut self, pos: usize, new_elem: T ) {
        self.vec.insert( pos, new_elem );
        self.size += 1;
    }

    /// Remove the last element from the vector and return it
    #[inline]
    pub fn pop(&mut self) -> T {
        let result = self.vec.pop();
        self.size -= 1;
        result.unwrap()
    }
}

impl<T: std::cmp::PartialEq> Vector<T> {
    /// Return the index of the first element in the Vector equal to a given value 
    #[inline]
    pub fn find(&self, value: T ) -> usize {
        let index = self.vec.iter().position( |x| *x == value );
        match index {
            None => panic!( "Entry not found in Vector." ),
            Some(index) => return index,
        }
    }
}

impl<T: Clone> Vector<T> {
    /// Create a new vector of specified size
    #[inline]
    pub fn new( size: usize, elem: T ) -> Self {
        let vec = vec![ elem; size ];
        Vector { vec, size }
    }

    /// Assign a value to every element in the vector
    #[inline]
    pub fn assign(&mut self, elem: T ) {
        for i in 0..self.size {
            self.vec[i] = elem.clone();
        }
    }
}

impl<T: Clone + Number> Vector<T> {
    /// Create a vector of zeros 
    #[inline]
    pub fn zeros( size: usize ) -> Self {
        let vec = vec![ T::zero(); size ];
        Vector{ vec, size }
    }

    /// Create a vector of ones
    #[inline]
    pub fn ones( size: usize ) -> Self {
        let vec = vec![ T::one(); size ];
        Vector{ vec, size }
    }
}

impl<T: Clone> Clone for Vector<T> {
    /// Clone the vector
    #[inline]
    fn clone(&self) -> Self {
        Self::create( self.vec.clone() )
    }
}

impl<T> Index<usize> for Vector<T> {
    type Output = T;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: usize ) -> &'a T {
        &self.vec[ index ]
    }
}

impl<T> IndexMut<usize> for Vector<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: usize ) -> &mut T {
        &mut self.vec[ index ] 
    }
}

impl<T: Clone + Neg<Output = T>> Neg for Vector<T> {
    type Output = Self;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        let mut result = self.vec.clone();
        for i in 0..result.len() {
            result[i] = -result[i].clone();
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> Add<Vector<T>> for Vector<T> {
    type Output = Self;
    /// Add the elements of two vectors together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        if self.size != plus.size { panic!( "Vector sizes do not agree (+)." ); }
        let mut result = Vec::new();
        for i in 0..self.size {
            result.push( self.vec[i].clone() + plus.vec[i].clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> Sub<Vector<T>> for Vector<T> {
    type Output = Self;
    /// Subtract the elements of one vector from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        if self.size != minus.size { panic!( "Vector sizes do not agree (-)." ); }
        let mut result = Vec::new();
        for i in 0..self.size {
            result.push( self.vec[i].clone() - minus.vec[i].clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> Mul<T> for Vector<T> {
    type Output = Self;
    /// Multiply a vector by a scalar (vector * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..self.size {
            result.push( self.vec[i].clone() * scalar.clone() );
        }
        Self::Output::create( result )
    }
}

impl Mul<Vector<f64>> for f64 {
    type Output = Vector<f64>;
    /// Allow multiplication on the left by f64 (f64 * vector)
    #[inline]
    fn mul(self, vector: Vector<f64>) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..vector.size {
            result.push( self.clone() * vector.vec[i].clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> Div<T> for Vector<T> {
    type Output = Self;
    /// Divide a vector by a scalar (vector / scalar)
    fn div(self, scalar: T) -> Self::Output {
        let mut result = Vec::new();
        for i in 0..self.size {
            result.push( self.vec[i].clone() / scalar.clone() );
        }
        Self::Output::create( result )
    }
}

impl<T: Clone + Number> AddAssign for Vector<T> {
    /// Add a vector to a mutable vector and assign the result ( += )
    fn add_assign(&mut self, rhs: Self) {
        if self.size != rhs.size { panic!( "Vector sizes do not agree (+=)." ); }
        for i in 0..self.size {
            self.vec[i] += rhs.vec[i].clone();
        }
    }
}

impl<T: Clone + Number> SubAssign for Vector<T> {
    /// Substract a vector from a mutable vector and assign the result ( -= )
    fn sub_assign(&mut self, rhs: Self) {
        if self.size != rhs.size { panic!( "Vector sizes do not agree (-=)." ); }
        for i in 0..self.size {
            self.vec[i] -= rhs.vec[i].clone();
        }
    }
}

impl<T: Clone + Number> AddAssign<T> for Vector<T> {
    /// Add the same value to every element in a mutable vector 
    fn add_assign(&mut self, rhs: T) {
        for i in 0..self.size {
            self.vec[i] += rhs.clone();
        }
    }
} 

impl<T: Clone + Number> SubAssign<T> for Vector<T> {
    /// Subtract the same value from every element in a mutable vector 
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..self.size {
            self.vec[i] -= rhs.clone();
        }
    }
} 

impl<T: Clone + Number> MulAssign<T> for Vector<T> {
    /// Multiply every element in a mutable vector by a scalar 
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.size {
            self.vec[i] *= rhs.clone();
        }
    }
} 

impl<T: Clone + Number> DivAssign<T> for Vector<T> {
    /// Divide every element in a mutable vector by a scalar 
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.size {
            self.vec[i] /= rhs.clone();
        }
    }
} 

impl<T: std::default::Default> Vector<T> {
    /// Resize the vector 
    #[inline]
    pub fn resize(&mut self, new_size: usize ) {
        self.vec.resize_with( new_size, Default::default );
        self.size = new_size;
    }
}

impl<T: Clone + Number> Vector<T> {
    /// Return the dot product of two vectors v.dot(w)
    #[inline]
    pub fn dot(&self, w: Vector<T>) -> T {
        if self.size != w.size { panic!( "Vector sizes do not agree dot()." ); }
        let mut result: T = T::zero();
        for i in 0..self.size {
            result += self.vec[i].clone() * w.vec[i].clone();
        }
        result
    }

    /// Return the sum of all the elements in the vector
    #[inline]
    pub fn sum(&self) -> T {
        self.sum_slice( 0, self.size - 1 )
    }

    /// Return the sum of the elements in the vector from index start to end 
    #[inline]
    pub fn sum_slice(&self, start: usize, end: usize) -> T {
        if start > end { panic!( "Vector sum: start > end." ); }
        if self.size <= start { panic!( "Vector range error." ); }
        if self.size <= end { panic!( "Vector range error." ); }
        let mut result: T = T::zero();
        for i in start..=end {
            result += self.vec[i].clone();
        }
        result
    }

    /// Return the product of all the elements in the vector
    #[inline]
    pub fn product(&self) -> T {
        self.product_slice( 0, self.size - 1 )
    }

    /// Return the product of the elements in the vector from index start to end
    #[inline]
    pub fn product_slice(&self, start: usize, end: usize) -> T {
        if start > end { panic!( "Vector sum: start > end." ); }
        if self.size <= start { panic!( "Vector range error." ); }
        if self.size <= end { panic!( "Vector range error." ); }
        let mut result: T = self.vec[start].clone();
        for i in start+1..=end {
            result *= self.vec[i].clone();
        }
        result
    }
}

impl<T: fmt::Display> Vector<T> {
    /// Print the vector to a file
    #[inline]
    pub fn output(&self, filename: &str) {
        let mut f = File::create(filename).expect("Unable to create file");
        for i in 0..self.size {                                                                                                                                                                  
            writeln!(f, "{}", self.vec[i]).unwrap();                                                                                                                            
        }
    }
}

impl<T: Clone + Signed> Vector<T> {
    /// Return a vector containing the absolute values of the elements 
    #[inline]
    pub fn abs(&self) -> Self {
        let size = self.size;
        let mut vec = vec![ T::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].abs();
        }
        Vector { vec, size }
    }

    /// Return the L1 norm: sum of absolute values 
    #[inline]
    pub fn norm_1(&self) -> T {
        let mut result = T::zero();
        for i in 0..self.size {
            result += self.vec[i].abs();
        }
        result
    }
}

impl<T: Clone + Signed> Vector<Complex::<T>> {
    /// Return a vector containing the complex conjugate of the elements
    #[inline]
    pub fn conj(&self) -> Self {
        let size = self.size;
        let mut vec = vec![ Complex::<T>::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().conj();
        }
        Vector { vec, size }
    }
}

impl<T: Clone + Number> Vector<Complex::<T>> {
    /// Return a vector containing the real parts of a vector of complex numbers
    #[inline]
    pub fn real(&self) -> Vector<T> {
        let size = self.size;
        let mut vec = vec![ T::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].clone().real;
        }
        Vector { vec, size }
    } 
}

impl Vector<f64> {
    /// Create a linearly spaced vector of f64 with n elements with lower limit a
    /// and upper limit b
    #[inline]
    pub fn linspace( a: f64, b: f64, size: usize ) -> Self {
        let mut vec = vec![ 0.0; size ];
        let h: f64 = ( b - a ) / ((size as f64) - 1.0);
        for i in 0..size {
            vec[i] = a + h * (i as f64);
        }
        Vector{ vec, size }
    }

    /// Create a nonuniform vector of f64 with n elements with lower limit a,
    /// upper limit b and exponent p (p=1 -> linear)
    #[inline]
    pub fn powspace( a: f64, b: f64, size: usize, p: f64 ) -> Self {
        let mut vec = vec![ 0.0; size ];
        for i in 0..size {
            vec[i] = a + (b - a) * libm::pow( (i as f64) / ((size as f64) - 1.0), p );
        }
        Vector{ vec, size }
    }

    /// Return the L2 norm: square root of the sum of the squares
    #[inline]
    pub fn norm_2(&self) -> f64 {
        let mut result = 0.0;
        for i in 0..self.size {
            result += libm::pow( self.vec[i].abs(), 2.0 );
        }
        libm::sqrt( result )
    }

    /// Return the Lp norm: p-th root of the sum of the absolute values 
    /// raised to the power p
    #[inline]
    pub fn norm_p(&self, p: f64 ) -> f64 {
        let mut result = 0.0;
        for i in 0..self.size {
            result += libm::pow( self.vec[i].abs(), p );
        }
        libm::pow( result, 1.0/p )
    }

    /// Return the Inf norm: largest absolute value element (p -> infinity)
    #[inline]
    pub fn norm_inf(&self) -> f64 {
        let mut result = self.vec[0].abs();
        for i in 1..self.size {
            if result < self.vec[i].abs() {
                result = self.vec[i].abs();
            }
        }
        result
    }

    // Create a vector containing n random elements between 0 and 1
    #[inline]
    pub fn random( size: usize ) -> Self {
        let mut vec = vec![ 0.0; size ];
        let mut rng = rand::thread_rng();
        for i in 0..size {
            vec[i] = rng.gen::<f64>()
        }
        Vector{ vec, size }
    }

    /// Return the dot product of two vectors v.dot(w)
    #[inline]
    pub fn dotf64(&self, w: &Vector<f64>) -> f64 {
        if self.size != w.size { panic!( "Vector sizes do not agree dot()." ); }
        let mut result = 0.0;
        for i in 0..self.size {
            result += self.vec[i] * w.vec[i];
        }
        result
    }
}

impl<T> fmt::Debug for Vector<T> where
    T: fmt::Debug
{
    /// Format the output [ v_0, v_1, ... ]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vec )
    }
}

impl<T> fmt::Display for Vector<T> where
    T: fmt::Debug
{
    /// Format the output [ v_0, v_1, ... ]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.vec )
    }
} 

impl<T: Ord + PartialEq + PartialOrd> Vector<T>
{
    /// Sort the Vector of elements
    pub fn sort(&mut self) {
        self.vec.sort_unstable();
    }
}

impl<T: PartialEq + PartialOrd> Vector<T>
{
    /// Sort the Vector of elements using a comparison function
    pub fn sort_by<F>(&mut self, compare: F )
    where
        F: FnMut(&T, &T) -> Ordering,
    {
        self.vec.sort_unstable_by( compare );
    }
}