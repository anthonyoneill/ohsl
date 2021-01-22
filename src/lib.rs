//! # Oak Hamilton Scientific Library
//!
//! `ohsl` is a collection of numerical routines and mathematical types
//! for use in scientific computing. 

//TODO quaternion.rs -> implement as a scalar and vector part in a struct + examples

pub mod elementary;
pub mod constant;
pub mod complex;
pub mod vector;
pub mod traits;
pub mod matrix;
pub mod mesh1d;
pub mod newton;
pub mod mesh2d;

// Re-exports
pub use self::complex::{Complex, Cmplx};
pub use self::vector::{Vector, Vec64};
pub use self::traits::{Number, Signed, Zero, One};
pub use self::matrix::{Matrix, Mat64, Sparse, Sparse64, Triplet, Tr64};
pub use self::mesh1d::Mesh1D;
pub use self::newton::Newton;
pub use self::mesh2d::Mesh2D;


// cargo test

#[cfg(test)]
mod tests {

    //pub use crate::complex::{Complex, Cmplx};
    //pub use crate::traits::{Zero, One};
    pub use crate::vector::{Vector, Vec64};
    pub use crate::matrix::{Matrix, Mat64, Sparse, Sparse64, Triplet, Tr64};
    //pub use crate::mesh1d::Mesh1D;
    //pub use crate::newton::Newton;
    //pub use crate::mesh2d::Mesh2D;

    #[test]
    fn test_example() {
        let v = Vector::<i32>::empty();
        assert_eq!( v.size(), 0 );
    }

}
 