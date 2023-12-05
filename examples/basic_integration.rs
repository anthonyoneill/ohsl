// cargo run --example basic_integration

extern crate ohsl;

pub use ohsl::vector::{Vector, Vec64};
pub use ohsl::mesh1d::Mesh1D;


fn main() {
    println!("--------- Basic integration ---------");

    println!( "  * Create a 1D mesh with 2 variables and set " );
    println!( "  * one equal 2x and the other to x^2. " );

    let nodes = Vec64::linspace( 0.0, 1.0, 101 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes.clone(), 2 );
    for i in 0..nodes.size() {
        let x = nodes[i];
        mesh[i][0] = 2.0 * x;
        mesh[i][1] = x * x;
    }

    println!( "  * number of nodes = {}", mesh.nnodes() );
    println!( "  * number of variables = {}", mesh.nvars() );
    println!( "  * Interpolate the value of each of the ");
    println!( "  * variables at x = 0.314 ");
    let vars = mesh.get_interpolated_vars( 0.314 );
    println!( "  * vars = {}", vars );

    println!( "  * Numerically integrate the variables over " );
    println!( "  * the domain (from 0 to 1) " );
    println!( "  * Integral 2x = {}", mesh.trapezium( 0 ) );
    println!( "  * Integral x^2 = {}", mesh.trapezium( 1 ) );

    // The mesh may be printed to a file using
    // mesh.output( "./output.txt", 5 );
    // here 5 is the precision of the output values.
    // Similarly a pre-existing file can be read into a mesh using
    // mesh.read( "./output.txt" );
    
    println!( "--- FINISHED ---" );
}

