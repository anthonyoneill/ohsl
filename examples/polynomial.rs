// cargo run --example polynomial

pub use ohsl::Polynomial;

fn main() {
    println!("----- Polynomials -----");

    // Create a new polynomial
    let p = Polynomial::<f64>::new( vec![ 1.0, 2.0, 3.0 ] );
    println!( "  * p(x) = {}", p );
    println!( "  * p[0] = {}", p[0] );
    println!( "  * p[1] = {}", p[1] );
    println!( "  * p[2] = {}", p[2] );
    println!( "  * p is of degree {}", p.degree().unwrap() );
    
    // Evaluate the polynomial at a given point
    println!( "  * p(0.5) = {}", p.eval( 0.5 ) );
    println!( "  * p(1.0) = {}", p.eval( 1.0 ) );

    // Determine the roots of the polynomial
    let roots = p.roots( true );
    println!( "  * p(x) = {} = 0 has {} roots", p, roots.size() );
    println!( "    * x0 = {:.6} {:.6} i", roots[0].real, roots[0].imag );
    println!( "    * x1 = {:.6} +{:.6} i", roots[1].real, roots[1].imag );
    
    // Arithmetic
    let q = Polynomial::cubic( 3.0, 5.0, 4.0, 0.0 );
    println!( "  * q(x) = {}", q );
    println!( "  * p(x) + q(x) = {}", &p + &q );
    println!( "  * p(x) - q(x) = {}", &p - &q );
    println!( "  * -p(x) = {}", -&p );
    println!( "  * p(x) * q(x) = {}", &p * &q );
    println!( "  * p(x) * 3  = {}", &p * 3.0 );
    let result = q.polydiv( &p );
    match result {
        Ok( ( q, r ) ) => println!( "  * q(x) / p(x) = {} remainder {}", q, r ),
        Err( e ) => println!( "  * q(x) / p(x) = {}", e )
    }
    
    // Derivatives
    println!( "  * p'(x) = {}", p.derivative() );
    println!( "  * p''(x) = {}", p.derivative_n( 2 ) );
    println!( "  * p'(0.5) = {}", p.derivative_at( 0.5, 1 ) );
    println!( "  * p''(0.5) = {}", p.derivative_at( 0.5, 2 ) );

    println!( "--- FINISHED ---" );
}