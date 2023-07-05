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

pub fn test_function<T>(input: T) -> T {
    input
}



// elementary.rs plan

// Could just use std f64 functions

/*
1. constants 
3. power + sqrt, cbrt
4. trigonometric functions
5. exponential and logarithmic functions
6. hyperbolic functions
7. Rounding, abs etc?
*/
