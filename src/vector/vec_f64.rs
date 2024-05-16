use rand::Rng;
pub use crate::vector::Vector;

impl Vector<f64> {
    /// Create a linearly spaced vector of f64 with n elements with lower limit a
    /// and upper limit b
    #[inline]
    pub fn linspace( a: f64, b: f64, size: usize ) -> Self {
        let mut vec = vec![ 0.0; size ];
        let h: f64 = ( b - a ) / ((size as f64) - 1.0);
        for i in 0..size {
            vec[i] = a + h * (i as f64);
        }
        Vector{ vec }
    }

    /// Create a nonuniform vector of f64 with n elements with lower limit a,
    /// upper limit b and exponent p (p=1 -> linear)
    #[inline]
    pub fn powspace( a: f64, b: f64, size: usize, p: f64 ) -> Self {
        let mut vec = vec![ 0.0; size ];
        for i in 0..size {
            vec[i] = a + (b - a) * f64::powf( (i as f64) / ((size as f64) - 1.0), p );
        }
        Vector{ vec }
    }

    /// Return the L2 norm: square root of the sum of the squares
    #[inline]
    pub fn norm_2(&self) -> f64 {
        let mut result = 0.0;
        for i in 0..self.size() {
            result += f64::powf( self.vec[i].abs(), 2.0 );
        }
        f64::sqrt( result )
    }

    /// Return the Lp norm: p-th root of the sum of the absolute values 
    /// raised to the power p
    #[inline]
    pub fn norm_p(&self, p: f64 ) -> f64 {
        let mut result = 0.0;
        for i in 0..self.size() {
            result += f64::powf( self.vec[i].abs(), p );
        }
        f64::powf( result, 1.0/p )
    }

    /// Return the Inf norm: largest absolute value element (p -> infinity)
    #[inline]
    pub fn norm_inf(&self) -> f64 {
        let mut result = self.vec[0].abs();
        for i in 1..self.size() {
            if result < self.vec[i].abs() {
                result = self.vec[i].abs();
            }
        }
        result
    }

    // Create a vector containing n random elements between 0 and 1
    #[inline]
    pub fn random( size: usize ) -> Self {
        let mut vec = vec![ 0.0; size ];
        let mut rng = rand::thread_rng();
        for i in 0..size {
            vec[i] = rng.gen::<f64>()
        }
        Vector{ vec }
    }
}
