pub use crate::vector::{Vector, Vec64};
pub use crate::mesh1d::Mesh1D;
//pub use crate::banded::Banded;
//pub use crate::matrix::Mat64;

pub use crate::traits::Number;

pub struct ODEIVP<T> {
    pub x0: f64,                    // Initial x value 
    pub xn: f64,                    // Final x value
    pub h: f64,                     // Step size (initial step size for adaptive methods)
    pub max_steps: usize,           // Maximum number of steps
    pub solution: Mesh1D<T, f64>,   // Solution mesh
    pub store_every: usize,         // Store every n steps for output
    pub tol: f64,                   // Tolerance for adaptive methods
}

impl<T: Clone + Number> ODEIVP<T> {
    /// Create a new ODEIVP solver 
    pub fn new( x0: f64, xn: f64, steps: usize ) -> Self {
        let h = (xn - x0) / (steps as f64);
        let max_steps: usize = steps;
        let solution = Mesh1D::<T, f64>::empty(); 
        ODEIVP { 
            x0, 
            xn, 
            h, 
            max_steps, 
            solution,
            store_every: 1,
            tol: 1.0e-8,
        }
    }

    /// Edit the maximum number of steps 
    pub fn max_steps(&mut self, max_steps: usize ) {
        self.max_steps = max_steps;
    }

    /// Edit the step size (initial step size for adaptive methods)
    pub fn step_size(&mut self, h: f64 ) {
        self.h = h;
    }

    /// Edit the store every n steps for output    #[inline]
    pub fn store_every(&mut self, store_every: usize ) {
        self.store_every = store_every;
    }

    /// Edit the tolerance for adaptive methods
    pub fn tolerance(&mut self, tol: f64 ) {
        self.tol = tol;
    }
}

impl ODEIVP<f64> {
    pub fn rk4(&mut self, 
        eqn: &dyn Fn(&Vec64, &f64) -> Vec64, 
        initial_conditions: Vec64
    ) -> Result<(), String> {
        let y_dash_0 = eqn( &initial_conditions, &self.x0 );
        if initial_conditions.size() != y_dash_0.size() {
            let msg = format!( 
                "ODEIVP error: initial conditions size ({}) does not match number of variables in the equation ({}).", 
                initial_conditions.size(), 
                y_dash_0.size() 
            );
            return Err( msg );
        }
        
        let mut x = self.x0;
        let h = self.h;
        let h2 = h / 2.0;
        let h6 = h / 6.0;
        let h3 = h / 3.0;
        let order = initial_conditions.size();

        let mut u = initial_conditions.clone();
        let mut z: Vec64;
        let (mut k1, mut k2, mut k3, mut k4): (Vec64, Vec64, Vec64, Vec64);

        let mut coords = Vec64::new( 1, x );
        let mut values: Vec<Vec64> = vec![ u.clone() ];

        for i in 0..self.max_steps {
            // k1 = F( u, x )
            k1 = eqn( &u, &x );
            z = &u + &(h2 * &k1);

            x += h2;

            // k2 = F( z, x + h/2 )
            k2 = eqn( &z, &x );
            z = &u + &(h2 * &k2);

            // k3 = F( z, x + h/2 )
            k3 = eqn( &z, &x );
            z = &u + &(h * &k3);

            x += h2;

            // k4 = F( z, x + h )
            k4 = eqn( &z, &x );
            u += h6 * &k1 + h3 * &k2 + h3 * &k3 + h6 * &k4;

            if i % self.store_every == 0 {
                coords.push( x );
                values.push( u.clone() );
            }

        }

        let size = coords.size();
        self.solution = Mesh1D::new( coords, order );
        for i in 0..size {
            self.solution.set_nodes_vars( i, values[i].clone() );
        }

        Ok(())
    }

    pub fn rkf45(&mut self, 
        eqn: &dyn Fn(&Vec64, &f64) -> Vec64, 
        initial_conditions: Vec64,
        //tol: f64,
        //TODO h_init?
    ) -> Result<(), String> {
        let y_dash_0 = eqn( &initial_conditions, &self.x0 );
        if initial_conditions.size() != y_dash_0.size() {
            let msg = format!( 
                "ODEIVP error: initial conditions size ({}) does not match number of variables in the equation ({}).", 
                initial_conditions.size(), 
                y_dash_0.size() 
            );
            return Err( msg );
        }

        let mut ok: bool;
        let mut step = 0;
        let mut x = self.x0;
        let mut h = self.h;
        let (mut c, mut diff): (f64, f64);

        const X2: f64 = 1. / 4.;
        const X3: f64 = 3. / 8.;
        const X4: f64 = 12. / 13.;
        const X5: f64 = 1.;
        const X6: f64 = 1. / 2.;

        const W21: f64 = 1. / 4.;

        const W31: f64 = 3. / 32.;
        const W32: f64 = 9. / 32.;

        const W41: f64 = 1932. / 2197.;
        const W42: f64 = -7200. / 2197.;
        const W43: f64 = 7296. / 2197.;

        const W51: f64 = 439. / 216.;
        const W52: f64 = -8.;
        const W53: f64 = 3680. / 513.;
        const W54: f64 = -845. / 4104.;

        const W61: f64 = -8. / 27.;
        const W62: f64 = 2.;
        const W63: f64 = -3544. / 2565.;
        const W64: f64 = 1859. / 4104.;
        const W65: f64 = -11. / 40.;

        const U1: f64 = 25. / 216.;
        const U3: f64 = 1408. / 2565.;
        const U4: f64 = 2197. / 4104.;
        const U5: f64 = -1. / 5.;

        const Z1: f64 = 16. / 135.;
        const Z3: f64 = 6656. / 12825.;
        const Z4: f64 = 28561. / 56430.;
        const Z5: f64 = -9. / 50.;
        const Z6: f64 = 2. / 55.;

        let order = initial_conditions.size();

        let mut u = initial_conditions.clone();
        let (mut z, mut e): (Vec64, Vec64);
        let (mut k1, mut k2, mut k3, mut k4, mut k5, mut k6): (Vec64, Vec64, Vec64, Vec64, Vec64, Vec64);

        let mut coords = Vec64::new( 1, x );
        let mut values: Vec<Vec64> = vec![ u.clone() ];

        loop {
            step += 1;

            // k1 = F( u, x )
            k1 = eqn( &u, &x );
            k1 *= h;
            z = &u + &(W21 * &k1);

            // k2 = F( z, x + X2*h )
            k2 = eqn( &z, &( x + X2 * h ) );
            k2 *= h;
            z = &u + &(W31 * &k1 + W32 * &k2);

            // k3 = F( z, x + X3*h )
            k3 = eqn( &z, &( x + X3 * h ) );
            k3 *= h;
            z = &u + &(W41 * &k1 + W42 * &k2 + W43 * &k3);

            // k4 = F( z, x + X4*h )
            k4 = eqn( &z, &( x + X4 * h ) );
            k4 *= h;
            z = &u + &(W51 * &k1 + W52 * &k2 + W53 * &k3 + W54 * &k4);

            // k5 = F( z, x + X5*h )
            k5 = eqn( &z, &( x + X5 * h ) );
            k5 *= h;
            z = &u + &(W61 * &k1 + W62 * &k2 + W63 * &k3 + W64 * &k4 + W65 * &k5);

            // k6 = F( z, x + X6*h )
            k6 = eqn( &z, &( x + X6 * h ) );
            k6 *= h;

            e = U1 * &k1 + U3 * &k3 + U4 * &k4 + U5 * &k5;
            z = Z1 * &k1 + Z3 * &k3 + Z4 * &k4 + Z5 * &k5 + Z6 * &k6;

            e -= &z;
            diff = e.norm_inf();
            c = (self.tol * h / ( 2. * diff )).sqrt().sqrt();
            ok = true;

            // is the first step ok? or does it need reducing?
            if (step == 1) && (c < 1.0) {
                // step needs reducing so start from initial value again
                ok = false;
                step = 1;
            }

            if ok {
                x += h;
                u += z;
                if step % self.store_every == 0 {
                    coords.push( x );
                    values.push( u.clone() );
                }
            }

            h *= c;

            if x + h > self.xn {
                h = self.xn - x;
            }

            if step >= self.max_steps {
                return Err( 
                    format!( "ODEIVP error: maximum number of steps ({}) exceeded.", self.max_steps ) 
                );
            }

            if !( (x - self.xn).abs() > self.tol ) {
                break;
            } 
        }

        let size = coords.size();
        self.solution = Mesh1D::new( coords, order );
        for i in 0..size {
            self.solution.set_nodes_vars( i, values[i].clone() );
        }

        Ok(())
    }
}