// cargo run --example sparse_laplace

extern crate ohsl;

pub use ohsl::vector::{Vector, Vec64};
pub use ohsl::matrix::{Matrix, Mat64};
use ohsl::sparse::Sparse;

use std::time::Instant;

fn main() {
    println!( "------------------ Laplace's equation ------------------" );
    println!( "  * Solve Laplace's equation on the unit square subject");
    println!( "  * to the given boundary conditions. The exact solution");
    println!( "  * is u(x,y) = y/[(1+x)^2 + y^2].");
    let start = Instant::now();
    let n: usize = 64;                          // Number of intervals
    let dx: f64 = 1.0 / ( n as f64 );           // Step size
    let size: usize = ( n + 1 ) * ( n + 1 );    // Size of the linear system
    let mut rhs = Vec64::zeros( size );         // Right hand side vector 

    let mut triplets = Vec::new();
    let mut row: usize = 0;
    
    // x = 0 boundary ( u = y / ( 1 + y^2 ) )
    let i = 0;
    for j in 0..n+1 {
        let y = (j as f64) * dx;
        triplets.push( ( row, i * ( n + 1 ) + j, 1.0 ) );
        rhs[ row ] = y / ( 1.0 + y * y );
        row += 1;
    }

    for i in 1..n {
        let x = (i as f64) * dx;
        // y = 0 boundary ( u = 0 )
        let j = 0;
        triplets.push( ( row, i * ( n + 1 ) + j, 1.0 ) );
        rhs[ row ] = 0.0;
        row += 1;
        // Interior nodes
        for j in 1..n {
            triplets.push( ( row, ( i - 1 ) * ( n + 1 ) + j, 1.0 ) );
            triplets.push( ( row, ( i + 1 ) * ( n + 1 ) + j, 1.0 ) );
            triplets.push( ( row, i * ( n + 1 ) + j - 1, 1.0 ) );
            triplets.push( ( row, i * ( n + 1 ) + j + 1, 1.0 ) );
            triplets.push( ( row, i * ( n + 1 ) + j, - 4.0 ) );
            rhs[ row ] = 0.0;
            row += 1;
        }
        // y = 1 boundary ( u = 1 / ( ( 1 + x )^2 + 1 ) )
        let j = n;
        triplets.push( ( row, i * ( n + 1 ) + j, 1.0 ) );
        rhs[ row ] = 1. / ( ( 1. + x ) * ( 1. + x ) + 1. );
        row += 1;
    }

    // x = 1 boundary ( u = y / ( 4 + y^2 ) )
    let i = n;
    for j in 0..n+1 {
        let y = (j as f64) * dx;
        triplets.push( ( row, i * ( n + 1 ) + j, 1.0 ) );
        rhs[ row ] = y / ( 4. + y * y );
        row += 1;
    }

    // Exact solution u = y / ( ( 1 + x )^2 + y^2 )
    let mut u_exact = Vec64::zeros( size );    
    row = 0;
    for i in 0..n+1 {
        let x = (i as f64) * dx;
        for j in 0..n+1 {
            let y = (j as f64) * dx;
            u_exact[ row ] = y / ( ( 1. + x ) * ( 1. + x ) + y * y );
            row += 1;
        }
    }
    // Create the sparse matrix from the triplets
    let sparse = Sparse::<f64>::from_triplets( size, size, &mut triplets );
    let duration = start.elapsed();
    println!("  * Time to create the matrix: {:?}", duration);
    // Solve the sparse system
    let mut u_guess = Vector::<f64>::random( size ); // A better guess would improve convergence
    let max_iter = 1000;
    let tol = 1e-8;
    let result = sparse.solve_qmr( &rhs, &mut u_guess, max_iter, tol );
    // Output time and error
    let duration = start.elapsed();
    match result {
        Ok( iter ) => println!( "  * The sparse system converged in {} iterations.", iter ),
        Err( error ) => println!( "  * The sparse system failed to converge, error = {}", error )
    }
    println!("  * Time to solve the system: {:?}", duration);
    let u_diff = u_guess - u_exact;
    println!("  * Solution error = {}", u_diff.norm_2() );
    println!( "-------------------- FINISHED ---------------------" );
}