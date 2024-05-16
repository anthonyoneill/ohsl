pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};
pub use crate::complex::Cmplx;

impl Matrix<f64> {
    /// Return the matrix one-norm (max absolute column sum)
    #[inline]
    pub fn norm_1(&self) -> f64 {
        let mut result: f64 = 0.0;
        for j in 0..self.cols {
            let mut sum: f64 = 0.0;
            for i in 0..self.rows {
                sum += self[(i,j)].abs();
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
                sum += self[(i,j)].abs();
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
                sum += f64::powf( self[(i,j)].abs(), p );
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
                result = result.max( self[(i,j)].abs() );
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

impl Matrix<Cmplx> {
    /// Create the Jacobian matrix of a complex vector valued function at a point
    /// using finite-differences
    #[inline]
    pub fn jacobian_cmplx( point: Vector<Cmplx>, func: &dyn Fn(Vector<Cmplx>) -> Vector<Cmplx>, delta: f64 ) -> Self {
        let n = point.size();
        let f = func( point.clone() );
        let m = f.size();
        let mut state = point.clone();
        let mut jac = Matrix::<Cmplx>::new( m, n, Cmplx::new( 0.0, 0.0 ) );
        for i in 0..n {
            state[i] += Cmplx::new( delta, 0.0 );
            let f_new = func( state.clone() ); 
            state[i] -= Cmplx::new( delta, 0.0 );
            jac.set_col( i, ( f_new - f.clone() ) / Cmplx::new( delta, 0.0 ) );
        }
        jac
    }
}