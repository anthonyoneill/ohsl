/// pi
pub const PI: f64 = 3.14159265358979323846;
/// pi/2 
pub const PI_2: f64 = 1.57079632679489661923;
/// Square root of 2
pub const SQRT2: f64 = 1.41421356237309504880;
/// e
pub const E: f64 = 2.71828182845904523536;

// Factorial function n! = n x (n-1) x (n-2) x ... x 2 x 1, n <= 20. 
// For larger values of n use the gamma function n! = gamma(n+1)
/*fn factorial(n: u64) -> f64 {
    if n > 20 {
        panic!("n must be <= 20, got {}.", n);
    }
    let p: u64 = (1..=n).product();
    p as f64 
}*/
// THIS IS A NAIVE IMPLEMENTATION SHOULD DO BETTER FOR PUBLIC FUNCTION 


// Continued fractions (generalised) ????
// Generalised continued fractions can be calculated using a recurrence 
// relation (look at the wikipedia page + numerical recipes p206 +
// book by Lorentzen and Waadeland ) 
// Could combine power series and continued fractions to calculate
// functions in different domains of convergence. 
// Ideally we want to remove the dependency upon libm 

// TODO - sqrtf64 - copying libm sqrt function (separate to general 
// sqrt function)

pub fn test_function<T>(input: T) -> T {
    input
}



// elementary.rs plan

/*
1. constants 
3. power + sqrt, cbrt
4. trigonometric functions
5. exponential and logarithmic functions
6. hyperbolic functions
7. Rounding, abs etc?
*/
