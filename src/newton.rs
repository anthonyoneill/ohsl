pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Complex;
pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};

pub struct Newton<T> {
    tol: f64,
    delta: f64,
    max_iter: usize,
    guess: T,
}

impl<T> Newton<T> {
    /// Create a new newton object
    #[inline]
    pub const fn new( guess: T ) -> Self {
        let tol: f64 = 1.0e-8;
        let delta: f64 = 1.0e-8;
        let max_iter: usize = 20;
        Newton { tol, delta, max_iter, guess }
    }

    /// Edit the convergence tolerance 
    #[inline]
    pub fn tolerance(&mut self, tolerance: f64 ) {
        self.tol = tolerance;
    }

    /// Edit the finite difference derivative step 
    #[inline]
    pub fn delta(&mut self, delta: f64 ) {
        self.delta = delta;
    }

    /// Edit the maximum number of iterations 
    #[inline]
    pub fn iterations(&mut self, iterations: usize ) {
        self.max_iter = iterations;
    }

    /// Edit the initial guess 
    #[inline]
    pub fn guess(&mut self, guess: T ) {
        self.guess = guess;
    }

}

impl<T: Copy> Newton<T> {
    /// Return the current parameters as a tuple 
    #[inline]
    pub fn parameters(&self) -> ( f64, f64, usize, T ) {
        let parameters = ( self.tol, self.delta, self.max_iter, self.guess );
        parameters
    }
}

impl Newton<f64> {
    /// Solve the equation via Newton iteration 
    #[inline]
    pub fn solve(&self, func: &dyn Fn(f64) -> f64 ) -> Result<f64, f64> {
        let mut current: f64 = self.guess;
        for _ in 0..self.max_iter {
            let deriv = ( func( current + self.delta ) - 
                          func( current - self.delta ) ) / ( 2.0 * self.delta );
            let dx = func(current) / deriv;
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( current );
            }
        }
        Err( current ) 
    }
}

impl Newton<Vec64> {
    /// Solve the vector equation via Newton iteration 
    #[inline] 
    pub fn solve(&self, func: &dyn Fn(Vec64) -> Vec64) -> Result<Vec64, Vec64> {
        let mut current: Vec64 = self.guess.clone();
        for _ in 0..self.max_iter {
            let f: Vec64 = func( current.clone() );
            let max_residual = f.norm_inf();
            let mut j = Mat64::jacobian( current.clone(), func, self.delta );
            let dx: Vec64 = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( current )
            }
        }
        Err( current )
    }

    /// Solve the vector equation via Newton iteration using the exact Jacobian
    #[inline] 
    pub fn solve_jacobian(&self, func: &dyn Fn(Vec64) -> Vec64, 
                                 jac: &dyn Fn(Vec64) -> Mat64 ) -> Result<Vec64, Vec64> {
        let mut current: Vec64 = self.guess.clone();
        for _ in 0..self.max_iter {
            let f: Vec64 = func( current.clone() );
            let max_residual = f.norm_inf();
            let mut j: Mat64 = jac( current.clone() ); 
            let dx: Vec64 = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( current )
            }
        }
        Err( current )
    }

}