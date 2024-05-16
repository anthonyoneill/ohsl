use core::ops::{Index, IndexMut};
pub use crate::vector::Vector;

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

impl<T> Vector<T> {
    /// Remove all of the elements from the vector
    #[inline]
    pub fn clear(&mut self) {
        self.vec.clear();
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
    }

    /// Push a new element onto the front of the vector
    #[inline]
    pub fn push_front(&mut self, elem: T ) {
        self.vec.insert( 0, elem );
    }

    /// Insert a new element into the Vector at a specified position 
    #[inline]
    pub fn insert(&mut self, pos: usize, new_elem: T ) {
        self.vec.insert( pos, new_elem );
    }

    /// Remove the last element from the vector and return it
    #[inline]
    pub fn pop(&mut self) -> T {
        let result = self.vec.pop();
        result.unwrap()
    }
}
