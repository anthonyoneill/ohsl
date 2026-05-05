pub use crate::traits::{Number, Signed, Zero, One};
pub use crate::complex::Cmplx;
pub use crate::vector::{Vector, Vec64};
pub use crate::matrix::{Matrix, Mat64};

pub struct Newton<T> {
    pub tol: f64,
    pub delta: f64,
    pub max_iter: usize,
    pub guess: T,
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
    pub fn solve(&self, func: &dyn Fn(&f64) -> f64 ) -> Result<(f64,usize), f64> {
        let mut current: f64 = self.guess;
        for iter in 0..self.max_iter {
            let f_plus = func( &(current + self.delta) );
            let f_minus = func( &(current - self.delta) );
            let dx = 2.0 * self.delta * func(&current) / ( f_plus - f_minus );
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current ) 
    }

    /// Solve the equation via Newton iteration using the exact derivative
    #[inline]
    pub fn solve_derivative(&self, 
        func: &dyn Fn(&f64) -> f64, 
        derivative: &dyn Fn(&f64) -> f64 
    ) -> Result<(f64,usize), f64> {
        let mut current: f64 = self.guess;
        for iter in 0..self.max_iter {
            let dx = func(&current) / derivative( &current );
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );   
            }
        }
        Err( current )
    }
}

impl Newton<Cmplx> {
    /// Solve the equation via Newton iteration 
    #[inline]
    pub fn solve(&self, func: &dyn Fn(&Cmplx) -> Cmplx ) -> Result<(Cmplx, usize), Cmplx> {
        let mut current: Cmplx = self.guess;
        let delta = Cmplx::new(self.delta, 0.0);
        for iter in 0..self.max_iter {
            let f_plus = func( &(current + delta) );
            let f_minus = func( &(current - delta) );
            let dx = 2.0 * delta * func(&current) / ( f_plus - f_minus );
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current ) 
    }

    /// Solve the equation via Newton iteration using the exact derivative
    #[inline]
    pub fn solve_derivative(&self, 
        func: &dyn Fn(&Cmplx) -> Cmplx, 
        derivative: &dyn Fn(&Cmplx) -> Cmplx 
    ) -> Result<(Cmplx, usize), Cmplx> {
        let mut current: Cmplx = self.guess;
        for iter in 0..self.max_iter {
            let dx = func(&current) / derivative( &current );
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current )
    }
}

impl Newton<Vec64> {
    /// Solve the vector equation via Newton iteration 
    #[inline] 
    pub fn solve(&self, func: &dyn Fn(&Vec64) -> Vec64) -> Result<(Vec64,usize), Vec64> {
        let mut current: Vec64 = self.guess.clone();
        for iter in 0..self.max_iter {
            let f: Vec64 = func( &current );
            let max_residual = f.norm_inf();
            let mut j = Mat64::jacobian( &current, func, self.delta );
            let dx: Vec64 = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( (current, iter+1) )
            }
        }
        Err( current )
    }

    /// Solve the vector equation via Newton iteration using the exact Jacobian
    #[inline] 
    pub fn solve_jacobian(&self, 
        func: &dyn Fn(&Vec64) -> Vec64, 
        jac: &dyn Fn(&Vec64) -> Mat64 
    ) -> Result<(Vec64,usize), Vec64> {
        let mut current: Vec64 = self.guess.clone();
        for iter in 0..self.max_iter {
            let f: Vec64 = func( &current );
            let max_residual = f.norm_inf();
            let mut j: Mat64 = jac( &current ); 
            let dx: Vec64 = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( (current, iter+1) )
            }
        }
        Err( current )
    }
}

impl Newton<Vector<Cmplx>> {
    /// Solve the vector equation via Newton iteration 
    #[inline] 
    pub fn solve(&self, 
        func: &dyn Fn(&Vector<Cmplx>) -> Vector<Cmplx>
    ) -> Result<(Vector<Cmplx>,usize), Vector<Cmplx>> {
        let mut current: Vector<Cmplx> = self.guess.clone();
        for iter in 0..self.max_iter {
            let f: Vector<Cmplx> = func( &current );
            let max_residual = f.norm_inf();
            let mut j = Matrix::jacobian_cmplx( &current, func, self.delta );
            let dx: Vector<Cmplx> = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( (current, iter+1) )
            }
        }
        Err( current )
    }

    /// Solve the vector equation via Newton iteration using the exact Jacobian
    #[inline] 
    pub fn solve_jacobian(&self, 
        func: &dyn Fn(&Vector<Cmplx>) -> Vector<Cmplx>, 
        jac: &dyn Fn(&Vector<Cmplx>) -> Matrix<Cmplx> 
    ) -> Result<(Vector<Cmplx>,usize), Vector<Cmplx>> {
        let mut current: Vector<Cmplx> = self.guess.clone();
        for iter in 0..self.max_iter {
            let f: Vector<Cmplx> = func( &current );
            let max_residual = f.norm_inf();
            let mut j: Matrix<Cmplx> = jac( &current ); 
            let dx: Vector<Cmplx> = j.solve_basic( &f );
            current -= dx;
            if max_residual <= self.tol {
                return Ok( (current, iter+1) )
            }
        }
        Err( current )
    }
}