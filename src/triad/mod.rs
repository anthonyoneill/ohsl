pub mod operations;

use std::fmt;

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

impl Triad<f64> {
    /// Create the Hessian of a vector valued function at a point using finite-differences 
    #[inline]
    pub fn hessian( point: &Vec64, func: &dyn Fn(&Vec64) -> Vec64, delta: f64 ) -> Self {
        let n = point.size();
        let f = func( point );
        let m = f.size();
        //TODO check m = n 
        let mut state = point.clone();
        let mut hes = Tri64::new( n, n, n, 0.0 );
        for j in 0..n {
            for k in 0..n {
                state[j] += delta;
                state[k] += delta;
                let f_pp = func( &state );
                state[k] -= 2. * delta;
                let f_pm = func( &state );
                state[j] -= 2. * delta;
                let f_mm = func( &state );
                state[k] += 2. * delta;
                let f_mp = func( &state );
                state[j] += delta;
                state[k] -= delta;
                for i in 0..m {
                    hes[( i, j, k )] = ( f_pp[i] - f_pm[i] - f_mp[i] + f_mm[i] ) / ( 4.0 * delta * delta );
                }
            }
        }
        hes
    }
}

impl<T> fmt::Debug for Triad<T> where
    T: fmt::Debug
{
    /// Format the output 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.panels {
            writeln!(f, "Panel {}:", i).unwrap();
            for j in 0..self.rows {
                for k in 0..self.cols {
                    write!(f, "\t{:?}", self[( i, j, k )] ).unwrap();
                }
                if j < self.rows-1 {
                    writeln!(f, "").unwrap();
                }
            }
            if i < self.panels-1 {
                writeln!(f, "").unwrap();
            }
        }
        write!(f, "")
    }
}