pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};

/* ----- Methods for solving linear systems ----- */

impl<T: Clone + Copy + Number + Signed + std::cmp::PartialOrd> Matrix<T> {
    #[inline]
    fn max_abs_in_column(&self, col: usize, start_row: usize) -> usize {
        let mut max_index: usize = 0;
        let mut max = T::zero();
        for i in start_row..self.rows {
            //if max < self[ i ][ col ].abs() {
            if max < self[(i,col)].abs() {
                //max = self[ i ][ col ].abs();
                max = self[(i,col)].abs();
                max_index = i;
            }
        }
        max_index
    }

    #[inline]
    fn backsolve(&self, x: &mut Vector<T> ) {
        let last = self.rows - 1;
        //x[ last ] = x[ last ]/ self[ last ][ last ];
        x[ last ] = x[ last ] / self[(last,last)];
        for n in 2..self.rows+1 {
            let k = self.rows - n;
            for j in self.rows-n+1..self.rows {
                let xj = x[ j ];
                //x[ k ] -= self[ k ][ j ] * xj;
                x[ k ] -= self[(k,j)] * xj;
            }
            //x[ k ] /= self[ k ][ k ];
            x[ k ] /= self[(k,k)];
        }
    }

    #[inline]
    fn partial_pivot(&mut self, x: &mut Vector<T>, k: usize ) {
        let pivot: usize = self.max_abs_in_column( k, k );
        self.swap_rows( pivot, k );
        x.swap( pivot, k );
    }

    #[inline]
    fn gauss_with_pivot(&mut self, x: &mut Vector<T> ){
        for k in 0..self.rows-1 {
            self.partial_pivot( x, k );
            for i in k+1..self.rows {
                //let elem = self[ i ][ k ] / self[ k ][ k ];
                let elem = self[(i,k)] / self[(k,k)];
                for j in k..self.rows {
                    //let kj = self[ k ][ j ];
                    let kj = self[(k,j)];
                    //self[ i ][ j ] -= elem * kj;
                    self[(i,j)] -= elem * kj;
                }
                let xk = x[ k ];
                x[ i ] -= elem * xk;
            }
        }
    }

    /// Solve the system of equations Ax=b where b is a specified vector 
    /// using Gaussian elimination (the matrix A is modified in the process)
    #[inline]
    pub fn solve_basic(&mut self, b: &Vector<T> ) -> Vector<T> {
        if self.rows != b.size() { panic!( "solve_basic error: rows != b.size()" ); }
        if self.rows != self.cols() { 
            panic!( "solve_basic error: matrix is not square" ); }
        let mut x: Vector<T> = b.clone();
        self.gauss_with_pivot( &mut x );
        self.backsolve( &mut x );
        x
    }

    /// Replace the matrix with its LU decomposition and return the number of pivots 
    /// and a permutation matrix
    #[inline]
    pub fn lu_decomp_in_place(&mut self) -> ( usize, Matrix<T> )  {
        if self.rows != self.cols() { panic!( "lu_decomp_in_place error: matrix is not square" ); }
        let mut pivots : usize = 0;
        let mut permutation = Matrix::<T>::eye( self.rows );
        for i in 0..self.rows() {
            let mut max_a = T::zero();
            let mut imax = i;
            for k in i..self.rows() {
                //let abs_a = self[ k ][ i ].abs();
                let abs_a = self[(k,i)].abs();
                if abs_a > max_a {
                    max_a = abs_a;
                    imax = k;
                }
            }
            //TODO check max_a to ensure matrix is not singular 
            if imax != i {
                permutation.swap_rows( i, imax );
                self.swap_rows( i, imax );
                pivots += 1;
            } 
            for j in i+1..self.rows() {
                //let ii = self[ i ][ i ];
                let ii = self[(i,i)];
                //self[ j ][ i ] /= ii;
                self[(j,i)] /= ii;
                for k in i+1..self.rows() { 
                    //let ji = self[ j ][ i ];
                    let ji = self[(j,i)];
                    //let ik = self[ i ][ k ];
                    let ik = self[(i,k)];
                    //self[ j ][ k ] -= ji * ik;
                    self[(j,k)] -= ji * ik;
                }
            }

        }
        ( pivots, permutation )
    }

    /// Solve the system of equations Ax=b where b is a specified vector 
    /// using LU decomposition (the matrix A is modified in the process)
    #[inline]
    pub fn solve_lu(&mut self, b: &Vector<T> ) -> Vector<T> {
        if self.rows != b.size() { panic!( "solve_LU error: rows != b.size()" ); }
        if self.rows != self.cols() { 
            panic!( "solve_LU error: matrix is not square" ); }
        let mut x: Vector<T> = b.clone();
        let ( _pivots, permutation ) = self.lu_decomp_in_place();
        x = permutation * x;
        for i in 0..self.rows() {
            for k in 0..i {
                let xk = x[ k ];
                //x[ i ] -= self[ i ][ k ] * xk;
                x[ i ] -= self[(i,k)] * xk;
            }
        }
        self.backsolve( &mut x );
        x
    }

    /// Calculate the determinant of the matrix ( via LU decomposition )
    #[inline]
    pub fn determinant(&self) -> T {
        let mut det = T::one();
        let mut temp = self.clone();
        let ( pivots, _permutation ) = temp.lu_decomp_in_place();
        for i in 0..self.rows() {
            //det *= temp.mat[ i ][ i ];
            det *= temp[(i,i)];
        }
        if pivots % 2 == 0 { det } else { - det }
    }

    /// Return the inverse of the matrix ( via LU decomposition )
    #[inline]
    pub fn inverse(&self) -> Matrix<T> {
        if self.rows != self.cols() { panic!( "inverse error: matrix is not square" ); }
        let mut lu: Matrix<T> = self.clone();
        let ( _pivots, mut inv ) = lu.lu_decomp_in_place();
        for j in 0..self.rows() {
            for i in 0..self.rows() {
                for k in 0..i {
                    //let inv_kj = inv.mat[ k ][ j ];
                    let inv_kj = inv[(k,j)];
                    //inv.mat[ i ][ j ] -=  lu.mat[ i ][ k ] * inv_kj;
                    inv[(i,j)] -= lu[(i,k)] * inv_kj;
                }
            }

            for i in (0..self.rows()).rev() {
                for k in i+1..self.rows() {
                    //let inv_kj = inv.mat[ k ][ j ];
                    let inv_kj = inv[(k,j)];
                    //inv.mat[ i ][ j ] -= lu.mat[ i ][ k ] * inv_kj;
                    inv[(i,j)] -= lu[(i,k)] * inv_kj;
                }
                //inv.mat[ i ][ j ] /= lu.mat[ i ][ i ];
                inv[(i,j)] /= lu[(i,i)];
            }
        }
        inv
    } 
}