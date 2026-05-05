pub mod operations;

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Complex;
pub use crate::vector::{Vector, Vec64};

#[derive(PartialEq)]
pub struct Triad<T> {
    tri: Vec<T>,
    panels: usize,
    rows: usize,
    cols: usize
}

pub type Tri64 = Triad<f64>;

impl<T> Triad<T> {
    /// Create a new triad of unspecified size
    #[inline]
    pub fn empty() -> Self {
        let tri = Vec::new();
        Triad { tri, panels: 0, rows: 0, cols: 0 }
    }

    /// Return the number of panels in the triad
    #[inline]
    pub fn panels(&self) -> usize {
        self.panels
    }

    /// Return the number of rows in the triad 
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Return the number of columns in the triad 
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Return the number of elements in the triad 
    #[inline]
    pub fn numel(&self) -> usize {
        self.cols * self.rows * self.panels
    }
}

impl<T: Clone + Number> Triad<T> {
    /// Create a new triad of specified size
    #[inline]
    pub fn new( panels: usize, rows: usize, cols: usize, elem: T ) -> Self {
        let size = panels * rows * cols;
        let mut tri = Vec::with_capacity( size );
        for _i in 0..size {
            tri.push( elem.clone() );
        }
        Triad { tri, panels, rows, cols }
    }

    /// Create a new triad from a vector of elements and specified size
    #[inline]
    pub fn from_vec( panels: usize, rows: usize, cols: usize, vec: Vec<T> ) -> Self {
        let size = panels * rows * cols;
        assert_eq!( vec.len(), size, 
            "The number of elements in the vector must match the specified size of the triad." 
        );
        Triad { tri: vec, panels, rows, cols }
    }
}

impl<T: Clone> Clone for Triad<T> {
    /// Clone the triad
    #[inline]
    fn clone(&self) -> Self {
        Triad { tri: self.tri.clone(), panels: self.panels, rows: self.rows, cols: self.cols }
    }
}