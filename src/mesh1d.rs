use core::ops::{Index, IndexMut};
use std::{fmt, fs::{File, read_to_string}, io::Write, str::FromStr}; 

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::vector::Vector;

pub struct Mesh1D<T, X> {
    nvars: usize,               // Number of variables
    nodes: Vector<X>,           // Vector for storing nodal points
    vars: Vec<Vector<T>>,       // Vec for storing variables a vector of variables at each node
}

impl<T: Clone + Number, X: Clone + Number + Copy> Mesh1D<T, X> {
    /// Create a new 1D mesh 
    #[inline]
    pub fn new( nodes: Vector<X>, nvars: usize ) -> Self {
        let node_vars = Vector::<T>::new( nvars, T::zero() );
        let mut vars = Vec::new();
        for _i in 0..nodes.size() {
            vars.push( node_vars.clone() );
        }
        Mesh1D { nvars, nodes, vars }
    }

    /// Return the number of nodal points in the mesh 
    #[inline]
    pub fn nnodes(&self) -> usize {
        self.nodes.size()
    }

    /// Return the number of variables stored at each node in the mesh 
    #[inline]
    pub fn nvars(&self) -> usize {
        self.nvars
    }

    /// Return the nodal coordinate at a specified index 
    #[inline]
    pub fn coord(&self, node: usize ) -> X {
        self.nodes[ node ]
    }

    /// Set the variables stored at a specified node
    #[inline]
    pub fn set_nodes_vars(&mut self, node: usize, vec: Vector<T> ) {
        //TODO node range checking
        if vec.size() != self.nvars { panic!( "Mesh1D error: set_nodes_vars " ); }
        self.vars[ node ] = vec;
    }

    /// Get the vector of variables stored at a specified node
    #[inline]
    pub fn get_nodes_vars(&self, node: usize ) -> Vector<T> {
        //TODO node range checking
        self.vars[ node ].clone()
    }

    /// Return the vector of nodal positions 
    #[inline]
    pub fn nodes(&self) -> Vector<X> {
        self.nodes.clone()
    }

}

impl Mesh1D<f64, f64> {
    /// Get the variables at an interpolated position ( first order scheme )
    #[inline]
    pub fn get_interpolated_vars(&self, x_pos: f64 ) -> Vector<f64> {
        //TODO range checking 
        let mut result = Vector::<f64>::new( self.nvars, 0.0 );
        for node in 0..self.nodes.size()-1 {
            if (( self.nodes[ node ] < x_pos ) && ( self.nodes[ node + 1 ] > x_pos ))
             || ( self.nodes[ node ] - x_pos ).abs() < 1.0e-7
             || ( self.nodes[ node + 1 ] - x_pos ).abs() < 1.0e-7
            {
                let delta_x: f64 = x_pos - self.nodes[ node ];
                let left = self.get_nodes_vars( node );
                let right = self.get_nodes_vars( node + 1 );
                let deriv = (right - left.clone()) / ( self.nodes[ node + 1 ] - self.nodes[ node ] );
                result = left + deriv * delta_x;
            }
        }
        result
    }

    /// Integrate a given variable over the domain (trapezium rule)
    #[inline]
    pub fn trapezium(&self, var: usize ) -> f64 {
        let mut sum: f64 = 0.0;
        for node in 0..self.nodes.size()-1 {
            let dx = self.nodes[ node + 1 ] - self.nodes[ node ];
            sum += 0.5 * dx * ( self.vars[ node ][ var ] 
                              + self.vars[ node + 1 ][ var ] );
        }
        sum
    }

    /// Read data from a file (overwrites nodes with file nodes)
    #[inline]
    pub fn read(&mut self, filename: &str ) {
        let data = read_to_string( filename ).expect( "Unable to read file");
        let split = data.split_whitespace();
        let vec: Vec<&str> = split.collect();

        let mut nodes = Vector::<f64>::empty();
        let node_vars = Vector::<f64>::new( self.nvars, 0.0 );

        for i in 0..vec.len() {
            if i % (self.nvars+1) == 0 {
                nodes.push( f64::from_str( vec[i] ).unwrap() );
            }  
        }
        self.vars.resize( nodes.size(), node_vars.clone() );
        self.nodes = nodes; 
        for i in 0..vec.len() {
            for var in 0..self.nvars {
                if i % (self.nvars+1) == var+1 {
                    self.vars[i / (self.nvars + 1)][var] = f64::from_str( vec[i] ).unwrap();
                }
            }
        }
    }

}

//TODO Mesh1D<Complex::<f64>, f64> get_interpolated_vars

impl<T, X> Index<usize> for Mesh1D<T, X> {
    type Output = Vector<T>;
    /// Indexing operator [] (read only) - returns the vector of variables
    /// stored at the specified node.
    #[inline]
    fn index<'a>(&'a self, node: usize ) -> &'a Vector<T> {
        &self.vars[ node ]
    }
}

impl<T, X> IndexMut<usize> for Mesh1D<T, X> {
    /// Indexing operator [] (read/write) - returns the vector of variables
    /// stored at the specified node.
    #[inline]
    fn index_mut(&mut self, node: usize ) -> &mut Vector<T> {
        &mut self.vars[ node ] 
    }
}


impl<T: fmt::Display, X: fmt::Display> Mesh1D<T, X> {
    /// Print the mesh to a file
    #[inline]
    pub fn output(&self, filename: &str, precision: usize ) {
        let mut f = File::create(filename).expect("Unable to create file");
        for i in 0..self.nodes.size() {  
            write!( f, "{number:.prec$} ", prec = precision, number = self.nodes[ i ] ).unwrap();
            for var in 0..self.nvars {
                write!( f, "{number:.prec$} ", prec = precision, number = self.vars[ i ][ var ] ).unwrap();
            }                                                                                                                                                               
            writeln!(f, "").unwrap();                                                                                                                           
        }
    }
} 