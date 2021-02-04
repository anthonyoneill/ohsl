// cargo run --example sparse_linear_system

extern crate ohsl;

pub use ohsl::vector::{Vector, Vec64};
pub use ohsl::matrix::{Matrix, Mat64, Sparse};

fn main() {
    println!( "------------------ Sparse linear system ------------------" );
    println!( "  * Solving the linear system Ax = b, for x, where" );
    let n: usize = 4;
    let mut a = Mat64::new( n, n, 0.0 );
    let mut triplets = Vec::new();
    for i in 0..n {
        a[i][i] = 2.0;
        triplets.push( ( i, i, 2.0 ) );
        if i > 0 {
            a[ i - 1 ][ i ] = 1.0;
            triplets.push( ( i - 1, i, 1.0 ) );
        }
        if i < n - 1 {
            a[ i + 1 ][ i ] = 1.0;
            triplets.push( ( i + 1, i, 1.0 ) );
        }
    }
    println!( "  * A ={}", a );
    let b = Vec64::random( n );
    println!( "  * b^T ={}", b );
    println!( "  * as both dense and sparse linear systems.");
    let dense_x = a.clone().solve_basic( b.clone() );
    println!( "  * The dense system gives the solution vector");
    println!( "  * x^T ={}", dense_x );
    let mut sparse = Sparse::from_triplets( triplets );
    let sparse_x = sparse.solve( b ).unwrap();
    println!( "  * The sparse system gives the solution vector");
    println!( "  * x^T ={}", sparse_x );
    let diff = dense_x - sparse_x;
    println!( "  * The maximum difference between the two is {:+.2e}", diff.norm_inf() );
    println!( "-------------------- FINISHED ---------------------" );
}