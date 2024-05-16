use std::cmp::Ordering;
pub use crate::traits::{Number, Signed};
pub use crate::vector::Vector;

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

impl<T: Copy> Vector<T> {
    /// Assign a value to every element in the vector
    #[inline]
    pub fn assign(&mut self, elem: T ) {
        for i in 0..self.size() {
            self.vec[i] = elem;
        }
    }
}

impl<T: std::default::Default> Vector<T> {
    /// Resize the vector 
    #[inline]
    pub fn resize(&mut self, new_size: usize ) {
        self.vec.resize_with( new_size, Default::default );
    }
}

impl<T: Copy + Number> Vector<T> {
    /// Return the dot product of two vectors v.dot(w)
    #[inline]
    pub fn dot(&self, w: &Vector<T>) -> T {
        if self.size() != w.size() { panic!( "Vector sizes do not agree dot()." ); }
        let mut result: T = T::zero();
        for i in 0..self.size() {
            result += self.vec[i] * w.vec[i];
        }
        result
    }

    /// Return the sum of all the elements in the vector
    #[inline]
    pub fn sum(&self) -> T {
        self.sum_slice( 0, self.size() - 1 )
    }

    /// Return the sum of the elements in the vector from index start to end 
    #[inline]
    pub fn sum_slice(&self, start: usize, end: usize) -> T {
        if start > end { panic!( "Vector sum: start > end." ); }
        if self.size() <= start { panic!( "Vector range error." ); }
        if self.size() <= end { panic!( "Vector range error." ); }
        let mut result: T = T::zero();
        for i in start..=end {
            result += self.vec[i].clone();
        }
        result
    }

    /// Return the product of all the elements in the vector
    #[inline]
    pub fn product(&self) -> T {
        self.product_slice( 0, self.size() - 1 )
    }

    /// Return the product of the elements in the vector from index start to end
    #[inline]
    pub fn product_slice(&self, start: usize, end: usize) -> T {
        if start > end { panic!( "Vector sum: start > end." ); }
        if self.size() <= start { panic!( "Vector range error." ); }
        if self.size() <= end { panic!( "Vector range error." ); }
        let mut result: T = self.vec[start].clone();
        for i in start+1..=end {
            result *= self.vec[i].clone();
        }
        result
    }
}

impl<T: Clone + Signed> Vector<T> {
    /// Return a vector containing the absolute values of the elements 
    #[inline]
    pub fn abs(&self) -> Self {
        let size = self.size();
        let mut vec = vec![ T::zero(); size ];
        for i in 0..size {
            vec[i] = self.vec[i].abs();
        }
        Vector { vec }//, size }
    }

    /// Return the L1 norm: sum of absolute values 
    #[inline]
    pub fn norm_1(&self) -> T {
        let mut result = T::zero();
        for i in 0..self.size() {
            result += self.vec[i].abs();
        }
        result
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
