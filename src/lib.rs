//! # Ostensibly Handy Scientific Library
//!
//! `ohsl` is a collection of numerical routines and mathematical types
//! for use in scientific computing. 

//TODO quaternion.rs -> implement as a scalar and vector part in a struct + examples

pub mod constant;
pub mod complex;
pub mod vector;
pub mod traits;
pub mod matrix;
pub mod mesh1d;
pub mod newton;
pub mod mesh2d;
pub mod sparse_matrix;
pub mod polynomial;
pub mod tridiagonal;
pub mod banded;

// Re-exports
pub use self::complex::{Complex, Cmplx};
pub use self::vector::{Vector, Vec64};
pub use self::traits::{Number, Signed, Zero, One};
pub use self::matrix::{Matrix, Mat64};
pub use self::mesh1d::Mesh1D;
pub use self::newton::Newton;
pub use self::mesh2d::Mesh2D;
pub use self::sparse_matrix::Sparse;
pub use self::polynomial::Polynomial;
pub use self::tridiagonal::Tridiagonal;
pub use self::banded::Banded;