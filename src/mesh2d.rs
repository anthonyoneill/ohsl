use core::ops::{Index, IndexMut};
use std::{fmt, fs::File, io::Write }; 

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::vector::Vector;
pub use crate::mesh1d::Mesh1D;
pub use crate::matrix::Matrix;

pub struct Mesh2D<T> {
    nvars: usize,               // Number of variables
    nx: usize,                  // Number of x_nodes
    ny: usize,                  // Number of y_nodes
    x_nodes: Vector<f64>,       // Vector for storing nodal points (x)
    y_nodes: Vector<f64>,       // Vector for storing nodal points (y)
    vars: Vec<Vector<T>>,       // Vec for storing variables a vector of variables at each node
}

impl<T: Clone + Number> Mesh2D<T> {
    /// Create a new 2D mesh 
    #[inline]
    pub fn new( x_nodes: Vector<f64>, y_nodes: Vector<f64>, nvars: usize ) -> Self {
        let node_vars = Vector::<T>::new( nvars, T::zero() );
        let mut vars = Vec::new();
        let nx = x_nodes.size();
        let ny = y_nodes.size();
        for _i in 0..nx {
            for _j in 0..ny {
                vars.push( node_vars.clone() );
            }
        }
        Mesh2D { nvars, nx, ny, x_nodes, y_nodes, vars }
    }

    /// Return the number of variables stored at each node in the mesh 
    #[inline]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// Return the number of nodal points in each direction 
    #[inline]
    pub fn nnodes(&self) -> ( usize, usize ) {
        ( self.nx, self.ny )
    }

    /// Return the spatial position of a specfied nodes as a tuple 
    #[inline]
    pub fn coord(&self, nodex: usize, nodey: usize ) -> ( f64, f64 ) {
        let px = self.x_nodes[ nodex ];
        let py = self.y_nodes[ nodey ];
        ( px, py )
    }

    /// Return the vector of x-nodal positions 
    #[inline]
    pub fn xnodes(&self) -> Vector<f64> {
        self.x_nodes.clone()
    }

    /// Return the vector of y-nodal positions 
    #[inline]
    pub fn ynodes(&self) -> Vector<f64> {
        self.y_nodes.clone()
    }

    /// Set the variables stored at a specified node 
    #[inline]
    pub fn set_nodes_vars(&mut self, nodex: usize, nodey: usize, vec: Vector<T> ) {
        if ( nodex > self.nx - 1 ) || ( nodey > self.ny - 1 ) { 
            panic!( "Mesh2D error: set_nodes_vars range error." ); 
        }
        if vec.size() != self.nvars { panic!( "Mesh2D error: set_nodes_vars nvars error." ); }
        self.vars[ nodex * self.ny + nodey ] = vec;
    }

    /// Get the vector of variables stored at a specified node
    #[inline]
    pub fn get_nodes_vars(&self, nodex: usize, nodey: usize ) -> Vector<T> {
        if ( nodex > self.nx - 1 ) || ( nodey > self.ny - 1 ) { 
            panic!( "Mesh2D error: get_nodes_vars range error." ); 
        }
        self.vars[ nodex * self.ny + nodey ].clone()
    }

    /// Assign an element to all entries in the mesh
    #[inline]
    pub fn assign(&mut self, element: T ) {
        for i in 0..self.nx {
            for j in 0..self.ny {
                for v in 0..self.nvars {
                    self.vars[ i * self.ny + j ][ v ] = element.clone();
                }
            }
        }
    }

    /// Return a cross section of the 2D mesh at a specified x node 
    #[inline]
    pub fn cross_section_xnode(&self, nodex: usize ) -> Mesh1D<T, f64> {
        let mut section = Mesh1D::<T, f64>::new( self.y_nodes.clone(), self.nvars );
        for nodey in 0..self.ny {
            section.set_nodes_vars( nodey, self.get_nodes_vars( nodex, nodey ) );
        }
        section
    }

    /// Return a cross section of the 2D mesh at a specified y node 
    #[inline]
    pub fn cross_section_ynode(&self, nodey: usize ) -> Mesh1D<T, f64> {
        let mut section = Mesh1D::<T, f64>::new( self.x_nodes.clone(), self.nvars );
        for nodex in 0..self.nx {
            section.set_nodes_vars( nodex, self.get_nodes_vars( nodex, nodey ) );
        }
        section
    }

    /// Return a matrix for a variable corresponding to each nodal point
    #[inline]
    pub fn var_as_matrix(&self, var: usize ) -> Matrix<T> {
        if var >= self.nvars { panic!( "Mesh2D var_as_matrix: index larger than # variables." ); }
        let mut m = Matrix::<T>::new( self.nx, self.ny, T::zero() );
        for i in 0..self.nx {
            for j in 0..self.ny {
                m[i][j] = self.vars[ i * self.ny + j ][ var ].clone();
            }
        }
        m
    } 

    /// Apply a function to the a specified variable in the mesh 
    #[inline]
    pub fn apply(&mut self, func: &dyn Fn(f64, f64) -> T, var: usize ) {
        for i in 0..self.nx {
            let x = self.x_nodes[i].clone();
            for j in 0..self.ny {
                let y = self.y_nodes[j].clone();
                self.vars[ i * self.ny + j ][ var ] = func( x, y );
            }
        }
    }

}

impl Mesh2D<f64> {
    /// Integrate a given variable over the domain (trapezium rule)
    #[inline]
    pub fn trapezium(&self, var: usize ) -> f64 {
        let mut sum: f64 = 0.0;
        for i in 0..self.nx-1 {
            let dx = self.x_nodes[ i + 1 ] - self.x_nodes[ i ];
            for j in 0..self.ny-1 {
                let dy = self.y_nodes[ j + 1 ] - self.y_nodes[ j ];
                sum += 0.25 * dx * dy * ( self.vars[ i * self.ny + j ][ var ]
                    + self.vars[ ( i + 1 ) * self.ny + j ][ var ]
                    + self.vars[ i * self.ny + j + 1 ][ var ]
                    + self.vars[ ( i + 1 ) * self.ny + j + 1 ][ var ] );
            }
        }
        sum
    }

    /// Integrate the square of a given variable over the domain (trapezium rule)
    #[inline]
    pub fn square_trapezium(&self, var: usize ) -> f64 {
        let mut sum: f64 = 0.0;
        for i in 0..self.nx-1 {
            let dx = self.x_nodes[ i + 1 ] - self.x_nodes[ i ];
            for j in 0..self.ny-1 {
                let dy = self.y_nodes[ j + 1 ] - self.y_nodes[ j ];
                sum += 0.25 * dx * dy * ( 
                      f64::powf( self.vars[ i * self.ny + j ][ var ].abs(), 2.0 )
                    + f64::powf( self.vars[ ( i + 1 ) * self.ny + j ][ var ].abs(), 2.0 )
                    + f64::powf( self.vars[ i * self.ny + j + 1 ][ var ].abs(), 2.0 )
                    + f64::powf( self.vars[ ( i + 1 ) * self.ny + j + 1 ][ var ].abs(), 2.0 ) );
            }
        }
        sum
    }
}

impl<T> Index<(usize, usize)> for Mesh2D<T> {
    type Output = Vector<T>;
    /// Indexing operator [] (read only) - returns the vector of variables
    /// stored at the node specified by a tuple.
    #[inline]
    fn index<'a>(&'a self, node: (usize, usize) ) -> &'a Vector<T> {
        &self.vars[ node.0 * self.ny + node.1 ]
    }
}

impl<T> IndexMut<(usize, usize)> for Mesh2D<T> {
    /// Indexing operator [] (read/write) - returns the vector of variables
    /// stored at the node specified by a tuple.
    #[inline]
    fn index_mut(&mut self, node: (usize, usize) ) -> &mut Vector<T> {
        &mut self.vars[ node.0 * self.ny + node.1 ] 
    }
}

impl<T: fmt::Display> Mesh2D<T> {
    /// Print the mesh to a file
    #[inline]
    pub fn output(&self, filename: &str, precision: usize ) {
        let mut f = File::create(filename).expect("Unable to create file");
        for j in 0..self.ny {
            for i in 0..self.nx {  
                write!( f, "{number:.prec$} ", prec = precision, number = self.x_nodes[ i ] ).unwrap();
                write!( f, "{number:.prec$} ", prec = precision, number = self.y_nodes[ j ] ).unwrap();
                for var in 0..self.nvars {
                    write!( f, "{number:.prec$} ", prec = precision, number = self.vars[ i * self.ny + j ][ var ] ).unwrap();
                }                                                                                                                                                               
                writeln!(f, "").unwrap();                                                                                                                           
            }
            writeln!(f, "").unwrap();
        }
    }

    /// Print the mesh of a single variable to a file 
    #[inline]
    pub fn output_var(&self, filename: &str, var: usize, precision: usize ) {
        let mut f = File::create(filename).expect("Unable to create file");
        for j in 0..self.ny {
            for i in 0..self.nx {  
                write!( f, "{number:.prec$} ", prec = precision, number = self.x_nodes[ i ] ).unwrap();
                write!( f, "{number:.prec$} ", prec = precision, number = self.y_nodes[ j ] ).unwrap();
                write!( f, "{number:.prec$} ", prec = precision, number = self.vars[ i * self.ny + j ][ var ] ).unwrap();                                                                                                                                                              
                writeln!(f, "").unwrap();                                                                                                                           
            }
            writeln!(f, "").unwrap();
        }
    }
    
} 
