pub use crate::complex::Cmplx;

pub struct Halley<T> {
    pub tol: f64,
    pub delta: f64,
    pub max_iter: usize,
    pub guess: T,
}

impl<T> Halley<T> {
    /// Create a new halley object
    #[inline]
    pub const fn new( guess: T ) -> Self {
        let tol: f64 = 1.0e-8;
        let delta: f64 = 1.0e-8;
        let max_iter: usize = 20;
        Halley { tol, delta, max_iter, guess }
    }
}

impl<T: Copy> Halley<T> {
    /// Return the current parameters as a tuple 
    #[inline]
    pub fn parameters(&self) -> ( f64, f64, usize, T ) {
        let parameters = ( self.tol, self.delta, self.max_iter, self.guess );
        parameters
    }
}

impl Halley<f64> {
    /// Solve the equation via Halley iteration 
    #[inline]
    pub fn solve(&self, func: &dyn Fn(&f64) -> f64 ) -> Result<(f64,usize), f64> {
        let mut current: f64 = self.guess;
        for iter in 0..self.max_iter {
            let f = func( &current );
            let f_plus = func( &(current + self.delta) );
            let f_minus = func( &(current - self.delta) );
            let df = f_plus - f_minus;
            let numerator = 2.0 * self.delta * f * df;
            let denominator = df * df - 2.0 * f * ( f_plus - 2.0 * f + f_minus );
            let dx = numerator / denominator;
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current ) 
    }

    /// Solve the equation via Halley iteration using the exact derivatives
    #[inline]
    pub fn solve_derivative(&self, 
        func: &dyn Fn(&f64) -> f64, 
        first: &dyn Fn(&f64) -> f64,
        second: &dyn Fn(&f64) -> f64
    ) -> Result<(f64,usize), f64> {
        let mut current: f64 = self.guess;
        for iter in 0..self.max_iter {
            let f = func( &current );
            let fd = first( &current );
            let fdd = second( &current );
            let numerator = f * fd; 
            let denominator = fd * fd - 0.5 * f * fdd;
            let dx = numerator / denominator;
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current )
    }
}

impl Halley<Cmplx> {
    /// Solve the equation via Halley iteration 
    #[inline]
    pub fn solve(&self, func: &dyn Fn(&Cmplx) -> Cmplx ) -> Result<(Cmplx,usize), Cmplx> {
        let mut current: Cmplx = self.guess;
        let delta = Cmplx::new(self.delta, 0.0);
        for iter in 0..self.max_iter {
            let f = func( &current );
            let f_plus = func( &(current + delta) );
            let f_minus = func( &(current - delta) );
            let df = f_plus - f_minus;
            let numerator = 2.0 * delta * f * df;
            let denominator = df * df - 2.0 * f * ( f_plus - 2.0 * f + f_minus );
            let dx = numerator / denominator;
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current ) 
    }

    /// Solve the equation via Halley iteration using the exact derivatives
    #[inline]
    pub fn solve_derivative(&self, 
        func: &dyn Fn(&Cmplx) -> Cmplx,
        first: &dyn Fn(&Cmplx) -> Cmplx,
        second: &dyn Fn(&Cmplx) -> Cmplx
    ) -> Result<(Cmplx,usize), Cmplx> {
        let mut current: Cmplx = self.guess;
        for iter in 0..self.max_iter {
            let f = func( &current );
            let fd = first( &current );
            let fdd = second( &current );
            let numerator = f * fd; 
            let denominator = fd * fd - 0.5 * f * fdd;
            let dx = numerator / denominator;
            current -= dx;
            if dx.abs() <= self.tol {
                return Ok( (current, iter+1) );
            }
        }
        Err( current )
    }
}