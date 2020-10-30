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
pub use self::matrix::{Matrix, Mat64};
pub use self::mesh1d::Mesh1D;
pub use self::newton::Newton;
pub use self::mesh2d::Mesh2D;


// cargo test

#[cfg(test)]
mod tests {

    //pub use crate::complex::{Complex, Cmplx};
    //pub use crate::traits::{Zero, One};
    pub use crate::vector::{Vector, Vec64};
    //pub use crate::matrix::{Matrix, Mat64};
    //pub use crate::mesh1d::Mesh1D;
    //pub use crate::newton::Newton;
    pub use crate::mesh2d::Mesh2D;

    #[test]
    fn test_example() {
        let v = Vector::<i32>::empty();
        assert_eq!( v.size(), 0 );
    }
    
    #[test]
    fn test_mesh2d_construction() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh[(0,0)][1] = 3.14;
        assert_eq!( mesh[(0,0)][1], 3.14 );
        assert_eq!( mesh.nvars(), 3 );
        assert_eq!( mesh.nnodes().0, 11 );
        assert_eq!( mesh.nnodes().1, 21 );
    }

    #[test]
    fn test_mesh2d_coord_nodes() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        assert_eq!( mesh.coord( 2, 4 ), ( 0.2, 0.4 ) );
        assert!( ( mesh.xnodes()[3] - 0.3 ).abs() < 1.0e-8 );
        assert!( ( mesh.ynodes()[19] - 1.9 ).abs() < 1.0e-8 );
    }

    #[test]
    fn test_mesh2d_set_nodes_vars() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        let vec = Vec64::new( 3, 3.14 );
        mesh.set_nodes_vars( 2, 4, vec );
        assert_eq!( mesh[(2,4)][0], 3.14 );
        assert_eq!( mesh[(2,4)][1], 3.14 );
        assert_eq!( mesh[(2,4)][2], 3.14 );
    }

    #[test]
    fn test_mesh2d_get_nodes_vars() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh[(0,0)][1] = 3.14;
        mesh[(0,0)][2] = 6.28;
        let vec = mesh.get_nodes_vars( 0, 0 );
        assert_eq!( vec[0], 0.0 );
        assert_eq!( vec[1], 3.14 );
        assert_eq!( vec[2], 6.28 );
    }

    #[test]
    fn test_mesh2d_assign() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh.assign( 2.718 );
        for i in 0..11 {
            for j in 0..21 {
                let vec = mesh.get_nodes_vars( i, j );
                assert_eq!( vec[0], 2.718 );
                assert_eq!( vec[1], 2.718 );
                assert_eq!( vec[2], 2.718 );
            }
        }
    }

    #[test]
    fn test_mesh2d_cross_section() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        let vec = Vec64::new( 3, 3.14 );
        mesh.set_nodes_vars( 2, 4, vec );
        let y_section = mesh.cross_section_xnode( 2 );
        assert_eq!( y_section.nnodes(), 21 );
        assert_eq!( y_section[4][0], 3.14 );
        assert_eq!( y_section[4][1], 3.14 );
        assert_eq!( y_section[4][2], 3.14 );
        let x_section = mesh.cross_section_ynode( 4 );
        assert_eq!( x_section.nnodes(), 11 );
        assert_eq!( x_section[2][0], 3.14 );
        assert_eq!( x_section[2][1], 3.14 );
        assert_eq!( x_section[2][2], 3.14 );
    }

    #[test]
    fn test_mesh2d_var_as_matrix() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 6 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 6 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh.assign( 2.718 );
        mesh[(0,0)][1] = 3.14;
        let mat = mesh.var_as_matrix( 1 );
        assert_eq!( mat[0][0], 3.14 );
        assert_eq!( mat[0][1], 2.718 );
    }

    /*#[test]
    fn test_mesh2d_output() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        for i in 0..mesh.xnodes().size() {
            let x = mesh.xnodes()[i].clone();
            for j in 0..mesh.ynodes().size() {
            let y = mesh.ynodes()[j].clone();
                mesh[(i,j)][0] = 2.0 * x * y;
                mesh[(i,j)][1] = x * x + y * y;
                mesh[(i,j)][2] = x - y;
            }
        }
        //mesh.output( "./output.txt", 5 );
        mesh.output_var( "./output.txt", 0, 5)
    }*/

    fn myfunction( x: f64, y: f64 ) -> f64 {
        2.0 * x * y
    }

    #[test]
    fn test_mesh2d_apply() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh.apply( &myfunction, 0 );
        assert!( ( mesh[(1,1)][0] - 0.02 ).abs() < 1.0e-8 );
    }

    #[test]
    fn test_mesh2d_trapezium() {
        let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
        let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
        let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
        mesh.apply( &myfunction, 0 );
        let integral = mesh.trapezium( 0 );
        assert!( ( integral - 2.0 ).abs() < 0.001 );
    }

}
 