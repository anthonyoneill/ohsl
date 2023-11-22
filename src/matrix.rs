use core::ops::{Index, IndexMut, Neg, Add, Sub, Mul, Div};
use core::ops::{AddAssign, SubAssign, MulAssign, DivAssign};
use std::{fmt, fs::File, io::Write, mem};//, cmp::Ordering};

pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Complex;
pub use crate::vector::{Vector, Vec64};

pub struct Matrix<T> {
    mat: Vec< Vector<T> >,
    rows: usize,
    cols: usize,
}

pub type Mat64 = Matrix<f64>;

impl<T> Matrix<T> {
    /// Create a new matrix of unspecified size
    #[inline]
    pub fn empty() -> Self {
        let row = Vector::<T>::empty();
        let mut mat = Vec::new();
        mat.push( row );
        let rows = 0;
        let cols = 0;
        Matrix { mat, rows, cols }
    }

    /// Create a matrix from an `std::vec::Vec<Vector<T>>`
    #[inline]
    pub fn create( mat: Vec< Vector<T> > ) -> Self {
        let rows = mat.len();
        let cols = mat[0].size();
        Matrix { mat, rows, cols }
    }

    /// Return the number of rows in the matrix 
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Return the number of columns in the matrix 
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Return the number of elements in the matrix 
    #[inline]
    pub fn numel(&self) -> usize {
        self.cols * self.rows
    }

    /// Remove all the elements from the matrix 
    #[inline]
    pub fn clear(&mut self) {
        self.mat.clear();
        self.rows = 0;
        self.cols = 0;
    }
}

impl<T: Clone + Number> Matrix<T> {
    /// Create a new matrix of specified size
    #[inline]
    pub fn new( rows: usize, cols: usize, elem: T ) -> Self {
        let row = Vector::<T>::new( cols, elem );
        let mut mat = Vec::new();
        for _i in 0..rows {
            mat.push( row.clone() );
        }
        Matrix { mat, rows, cols }
    }

    /// Get a row of the matrix as a vector
    #[inline]
    pub fn get_row(&self, row: usize ) -> Vector<T> {
        //if row < 0 || self.rows <= row { panic!( "Matrix range error in get_row" ); }
        self.mat[ row ].clone()
    }

    /// Get a column of the matrix as a vector 
    #[inline]
    pub fn get_col(&self, col: usize ) -> Vector<T> {
        //if col < 0 || self.rows <= col { panic!( "Matrix range error in get_col" ); }
        let mut result = Vector::<T>::new( self.rows, T::zero() );
        for i in 0..self.rows {
            result[ i ] = self[i][col].clone()
        }
        result
    }

    /// Set a row of the matrix using a vector 
    #[inline]
    pub fn set_row(&mut self, row: usize, vec: Vector<T> ) {
        //if vec.size() != self.cols { panic!( "Matrix size error in set_row" ); }
        //if row < 0 || self.rows <= row { panic!( "Matrix range error in set_row" ); }
        self[ row ] = vec;
    }

    /// Set a column of the matrix using a vector 
    #[inline]
    pub fn set_col(&mut self, col: usize, vec: Vector<T> ) {
        //if vec.size() != self.rows { panic!( "Matrix size error in set_col" ); }
        //if col < 0 || self.rows <= col { panic!( "Matrix range error in set_col" ); }
        for i in 0..self.rows {
            self[i][col] = vec[i].clone();
        }
    }

    /// Delete a row from the matrix 
    #[inline]
    pub fn delete_row(&mut self, row: usize ) {
        //if row < 0 || self.rows <= row { panic!( "Matrix range error in delete_row" ); }
        self.mat.remove( row );
        self.rows -= 1;
    }

    /// Multiply the matrix by a (column) vector and return a vector 
    #[inline]
    pub fn multiply(&self, vec: Vector<T> ) -> Vector<T> {
        if vec.size() != self.cols { panic!( "Matrix dimensions do not agree in multiply." ); }
        let mut result = Vector::<T>::empty();
        for row in 0..self.rows {
           result.push( self[row].dot( vec.clone() ) );
        }
        result
    }

    /// Create a square identity matrix of specified size 
    #[inline]
    pub fn eye( size: usize ) -> Self {
        let mut identity = Matrix::<T>::new( size, size, T::zero() ); 
        for i in 0..size {
            identity[i][i] = T::one();
        }
        identity
    }

    /// Resize the matrix (empty entries are appended if necessary)
    #[inline]
    pub fn resize(&mut self, n_rows: usize, n_cols: usize ) {
        let temp = self.clone();
        *self = Matrix::<T>::new( n_rows, n_cols, T::zero() );
        for i in 0..n_rows {
            for j in 0..n_cols {
                if i < temp.rows() && j < temp.cols() {
                    self[i][j] = temp[i][j].clone();
                }
            }
        }
    }

    /// Transpose the matrix in place 
    #[inline]
    pub fn transpose_in_place(&mut self) {
        if self.rows == self.cols {
            for i in 0..self.rows {
                for j in i+1..self.cols {
                    /*let temp = self[i][j].clone();
                    self[i][j] = self[j][i].clone();
                    self[j][i] = temp;*/

                    let mut temp = self[i][j].clone();
                    mem::swap( &mut self[j][i], &mut temp );
                    self[i][j] = temp;
                }
            }
        } else {
            let row = Vector::<T>::new( self.rows, T::zero() );
            let mut temp = Vec::new();
            for _i in 0..self.cols {
                temp.push( row.clone() );
            }
            for i in 0..self.rows {
                for j in 0..self.cols {
                    temp[j][i] = self[i][j].clone();
                }
            }
            self.mat = temp;
            mem::swap( &mut self.rows, &mut self.cols );
        }
    }

    /// Return the transpose of the matrix 
    #[inline]
    pub fn transpose(&self) -> Matrix<T> {
        let mut temp: Matrix<T> = self.clone();
        temp.transpose_in_place();
        temp
    }

    /// Swap two rows of the matrix 
    #[inline]
    pub fn swap_rows(&mut self, row_1: usize, row_2: usize ) {
        if self.rows <= row_1 || self.rows <= row_2 { 
            panic!( "Matrix swap row range error." ); 
        } 
        let mut temp = self.mat[ row_1 ].clone();
        mem::swap( &mut self.mat[ row_2 ], &mut temp );
        self.mat[ row_1 ] = temp;
    }

    /// Swap two elements of the matrix 
    #[inline]
    pub fn swap_elem(&mut self, row_1: usize, col_1: usize, row_2: usize, col_2: usize ) {
        let mut temp = self[row_1][col_1].clone();
        mem::swap( &mut self[row_2][col_2], &mut temp );
        self[row_1][col_1] = temp;
    }

    /// Fill the matrix with specified elements
    #[inline]
    pub fn fill(&mut self, elem: T ) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] = elem.clone();
            }
        }
    }

    /// Fill the leading diagonal of the matrix with specified elements
    #[inline]
    pub fn fill_diag(&mut self, elem: T ) {
        let n: usize = if self.cols < self.rows { self.cols } else { self.rows };
        for i in 0..n {
            self[i][i] = elem.clone();
        }
        
    }

    /// Fill a diagonal band of the matrix with specified elements
    /// ( offset above main diagonal +, below main diagonal - )
    #[inline]
    pub fn fill_band(&mut self, offset: isize, elem: T ) {
        for row in 0..self.rows {
            let i = (row as isize) + offset; 
            if (i as usize) < self.cols &&  i >= 0 {
                self[ row ][ i as usize ] = elem.clone();
            } 
        }
    }

    /// Fill the main three diagonals of the matrix with specified elements 
    #[inline]
    pub fn fill_tridiag(&mut self, lower: T, diag: T, upper: T ) {
        self.fill_band( -1, lower );
        self.fill_diag( diag );
        self.fill_band( 1, upper );
    }

}

impl<T: Clone + Copy + Number + Signed + std::cmp::PartialOrd> Matrix<T> {

    /* ----- Methods for solving linear systems ----- */

    #[inline]
    fn max_abs_in_column(&self, col: usize, start_row: usize) -> usize {
        let mut max_index: usize = 0;
        let mut max = T::zero();
        for i in start_row..self.rows {
            if max < self[ i ][ col ].abs() {
                max = self[ i ][ col ].abs();
                max_index = i;
            }
        }
        max_index
    }

    #[inline]
    fn backsolve(&self, x: &mut Vector<T> ) {
        let last = self.rows - 1;
        x[ last ] = x[ last ]/ self[ last ][ last ];
        for n in 2..self.rows+1 {
            let k = self.rows - n;
            for j in self.rows-n+1..self.rows {
                let xj = x[ j ];
                x[ k ] -= self[ k ][ j ] * xj;
            }
            x[ k ] /= self[ k ][ k ];
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
                let elem = self[ i ][ k ] / self[ k ][ k ];
                for j in k..self.rows {
                    let kj = self[ k ][ j ];
                    self[ i ][ j ] -= elem * kj;
                }
                let xk = x[ k ];
                x[ i ] -= elem * xk;
            }
        }
    }

    /// Solve the system of equations Ax=b where b is a specified vector 
    /// using Gaussian elimination (the matrix A is modified in the process)
    #[inline]
    pub fn solve_basic(&mut self, b: Vector<T> ) -> Vector<T> {
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
        if self.rows != self.cols() { 
            panic!( "lu_decomp_in_place error: matrix is not square" ); }
        let mut pivots : usize = 0;
        let mut permutation = Matrix::<T>::eye( self.rows );
        for i in 0..self.rows() {
            let mut max_a = T::zero();
            let mut imax = i;
            for k in i..self.rows() {
                let abs_a = self[ k ][ i ].abs();
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
                let ii = self[ i ][ i ];
                self[ j ][ i ] /= ii;
                for k in i+1..self.rows() { 
                    let ji = self[ j ][ i ];
                    let ik = self[ i ][ k ];
                    self[ j ][ k ] -= ji * ik;
                }
            }

        }
        ( pivots, permutation )
    }

    /// Solve the system of equations Ax=b where b is a specified vector 
    /// using LU decomposition (the matrix A is modified in the process)
    #[inline]
    pub fn solve_lu(&mut self, b: Vector<T> ) -> Vector<T> {
        if self.rows != b.size() { panic!( "solve_LU error: rows != b.size()" ); }
        if self.rows != self.cols() { 
            panic!( "solve_LU error: matrix is not square" ); }
        let mut x: Vector<T> = b.clone();
        let ( _pivots, permutation ) = self.lu_decomp_in_place();
        x = permutation * x;
        for i in 0..self.rows() {
            for k in 0..i {
                let xk = x[ k ];
                x[ i ] -= self[ i ][ k ] * xk;
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
            det *= temp.mat[ i ][ i ];
        }
        if pivots % 2 == 0 {
            det
        } else {
            - det
        }
    }

    /// Return the inverse of the matrix ( via LU decomposition )
    #[inline]
    pub fn inverse(&self) -> Matrix<T> {
        if self.rows != self.cols() { 
            panic!( "inverse error: matrix is not square" ); }
        let mut lu: Matrix<T> = self.clone();
        let ( _pivots, mut inv ) = lu.lu_decomp_in_place();
        for j in 0..self.rows() {
            for i in 0..self.rows() {
                for k in 0..i {
                    let inv_kj = inv.mat[ k ][ j ];
                    inv.mat[ i ][ j ] -=  lu.mat[ i ][ k ] * inv_kj;
                }
            }

            for i in (0..self.rows()).rev() {
                for k in i+1..self.rows() {
                    let inv_kj = inv.mat[ k ][ j ];
                    inv.mat[ i ][ j ] -= lu.mat[ i ][ k ] * inv_kj;
                }
                inv.mat[ i ][ j ] /= lu.mat[ i ][ i ];
            }
        }
        inv
    }
    
}

impl Matrix<f64> {
    /// Return the matrix one-norm (max absolute column sum)
    #[inline]
    pub fn norm_1(&self) -> f64 {
        let mut result: f64 = 0.0;
        for j in 0..self.cols {
            let mut sum: f64 = 0.0;
            for i in 0..self.rows {
                sum += self[i][j].abs()
            }
            result = result.max( sum )
        }
        result
    }

    /// Return the matrix inf-norm (max absolute row sum)
    #[inline]
    pub fn norm_inf(&self) -> f64 {
        let mut result: f64 = 0.0;
        for i in 0..self.rows {
            let mut sum: f64 = 0.0;
            for j in 0..self.cols {
                sum += self[i][j].abs()
            }
            result = result.max( sum )
        }
        result
    }

    /// Return the matrix p-norm (p=2 is Frobenius, p=inf is max norm)
    #[inline]
    pub fn norm_p(&self, p: f64 ) -> f64 {
        let mut sum: f64 = 0.0;
        for i in 0..self.rows {
            for j in 0..self.cols {
                sum += f64::powf( self[i][j].abs(), p );
            }
        }
        f64::powf( sum, 1.0/p )
    }

    /// Return the matrix Frobenius norm 
    #[inline]
    pub fn norm_frob(&self) -> f64 {
        self.norm_p( 2.0 )
    }

    /// Return the entrywise max-norm of the matrix 
    #[inline]
    pub fn norm_max(&self) -> f64 {
        let mut result: f64 = 0.0;
        for i in 0..self.rows {
            for j in 0..self.cols {
                result = result.max( self[i][j].abs() );
            }
        }
        result
    }

    /// Create the Jacobian matrix of a vector valued function at a point
    /// using finite-differences 
    #[inline]
    pub fn jacobian( point: Vec64, func: &dyn Fn(Vec64) -> Vec64, delta: f64 ) -> Self {
        let n = point.size();
        let f = func( point.clone() );
        let m = f.size();
        let mut state = point.clone();
        let mut jac = Mat64::new( m, n, 0.0 );
        for i in 0..n {
            state[i] += delta;
            let f_new = func( state.clone() ); 
            state[i] -= delta;
            jac.set_col( i, ( f_new - f.clone() ) / delta );
        }
        jac
    }

}

impl<T: Clone> Clone for Matrix<T> {
    /// Clone the matrix
    #[inline]
    fn clone(&self) -> Self {
        Self::create( self.mat.clone() )
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = Vector<T>;
    /// Indexing operator [] (read only)
    #[inline]
    fn index<'a>(&'a self, index: usize ) -> &'a Vector<T> {
        &self.mat[ index ]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    /// Indexing operator [] (read/write)
    #[inline]
    fn index_mut(&mut self, index: usize ) -> &mut Vector<T> {
        &mut self.mat[ index ] 
    }
}

impl<T: Clone + Neg<Output = T>> Neg for Matrix<T> {
    type Output = Self;
    /// Return the unary negation ( unary - ) of each element
    #[inline]
    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = -result[i][j].clone();
            }
        }
        result
    }
}

impl<T: Clone + Number> Add<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Add the elements of two matrices together ( binary + )
    #[inline]
    fn add(self, plus: Self) -> Self::Output {
        if self.rows != plus.rows { panic!( "Matrix row dimensions do not agree (+)." ); }
        if self.cols != plus.cols { panic!( "Matrix col dimensions do not agree (+)." ); }
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j].clone() + plus[i][j].clone();
            }
        }
        result
    }
}

impl<T: Clone + Number> Sub<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Subtract the elements of one matrix from another ( binary - )
    #[inline]
    fn sub(self, minus: Self) -> Self::Output {
        if self.rows != minus.rows { panic!( "Matrix row dimensions do not agree (-)." ); }
        if self.cols != minus.cols { panic!( "Matrix col dimensions do not agree (-)." ); }
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j].clone() - minus[i][j].clone();
            }
        }
        result
    }
}

impl<T: Clone + Number> Mul<T> for Matrix<T> {
    type Output = Self;
    /// Multiply a matrix by a scalar (matrix * scalar)
    #[inline]
    fn mul(self, scalar: T) -> Self::Output {
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j].clone() * scalar.clone();
            }
        }
        result
    }
}

impl Mul<Matrix<f64>> for f64 {
    type Output = Matrix<f64>;
    /// Allow multiplication on the left by f64 (f64 * matrix)
    #[inline]
    fn mul(self, matrix: Matrix<f64>) -> Self::Output {
        let mut result = Matrix::<f64>::new( matrix.rows(), matrix.cols(), 0.0 );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = matrix[i][j].clone() * self.clone();
            }
        }
        result
    }
}

impl<T: Clone + Number> Div<T> for Matrix<T> {
    type Output = Self;
    /// Divide a matrix by a scalar (matrix / scalar)
    fn div(self, scalar: T) -> Self::Output {
        let mut result = Matrix::<T>::new( self.rows(), self.cols(), T::zero() );
        for i in 0..result.rows() {
            for j in 0..result.cols() {
                result[i][j] = self[i][j].clone() / scalar.clone();
            }
        }
        result
    }
}

impl<T: Clone + Number> AddAssign for Matrix<T> {
    /// Add a matrix to a mutable matrix and assign the result ( += )
    fn add_assign(&mut self, rhs: Self) {
        if self.rows != rhs.rows { panic!( "Matrix row dimensions do not agree (+=)." ); }
        if self.cols != rhs.cols { panic!( "Matrix col dimensions do not agree (+=)." ); }
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] += rhs[i][j].clone();
            }
        }
    }
}

impl<T: Clone + Number> SubAssign for Matrix<T> {
    /// Subtract a matrix from a mutable matrix and assign the result ( -= )
    fn sub_assign(&mut self, rhs: Self) {
        if self.rows != rhs.rows { panic!( "Matrix row dimensions do not agree (-=)." ); }
        if self.cols != rhs.cols { panic!( "Matrix col dimensions do not agree (-=)." ); }
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] -= rhs[i][j].clone();
            }
        }
    }
}

impl<T: Clone + Number> MulAssign<T> for Matrix<T> {
    /// Multiply a mutable matrix by a scalar (matrix *= scalar)
    fn mul_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] *= rhs.clone();
            }
        }
    }
} 

impl<T: Clone + Number> DivAssign<T> for Matrix<T> {
    /// Divide a mutable matrix by a scalar (matrix /= scalar)
    fn div_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] /= rhs.clone();
            }
        }
    }
}

impl<T: Clone + Number> AddAssign<T> for Matrix<T> {
    /// Add the same value to every element in a mutable matrix
    fn add_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] += rhs.clone();
            }
        }
    }
}

impl<T: Clone + Number> SubAssign<T> for Matrix<T> {
    /// Subtract the same value from every element in a mutable matrix
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self[i][j] -= rhs.clone();
            }
        }
    }
}

impl<T: Clone + Number> Mul<Matrix<T>> for Matrix<T> {
    type Output = Self;
    /// Multiply two matrices together ( matrix * matrix )
    #[inline]
    fn mul(self, mul: Self) -> Self::Output {
        if self.cols != mul.rows { panic!( "Matrix dimensions do not agree (*)." ); }
        let mut result = Matrix::<T>::new( self.rows(), mul.cols(), T::zero() );
        for col in 0..mul.cols() {
            result.set_col( col, self.multiply( mul.get_col( col ) ) );
        }
        result
    }
}

impl<T: Clone + Number> Mul<Vector<T>> for Matrix<T> {
    type Output = Vector<T>;
    /// Multiply a matrix with a (column) vector ( matrix * vector )
    #[inline]
    fn mul(self, vec: Vector<T> ) -> Vector<T> {
        self.multiply( vec )
    }
}

impl<T> fmt::Debug for Matrix<T> where
    T: fmt::Debug
{
    /// Format the output 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows-1 {
            writeln!(f, "\t{:?}", self.mat[i] ).unwrap();
        }
        write!(f, "\t{:?}", self.mat[self.rows-1] )
    }
}

impl<T> fmt::Display for Matrix<T> where
    T: fmt::Debug
{
    /// Format the output 
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows-1 {
            writeln!(f, "\t{:?}", self.mat[i] ).unwrap();
        }
        write!(f, "\t{:?}", self.mat[self.rows-1] )
    }
} 

impl<T: fmt::Display> Matrix<T> {
    /// Print the matrix to a file
    #[inline]
    pub fn output(&self, filename: &str) {
        let mut f = File::create(filename).expect("Unable to create file");
        for i in 0..self.rows {  
            for j in 0..self.cols {
                write!(f, "\t{}", self[i][j] ).unwrap();
            }                                                                                                                                                                
            writeln!(f, "").unwrap();                                                                                                                            
        }
    }
}