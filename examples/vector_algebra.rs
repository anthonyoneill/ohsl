// cargo run --example vector_algebra

extern crate ohsl;

pub use ohsl::vector::{Vector, Vec64};
pub use ohsl::constant::*;


fn main() {
    println!("----- Vector algebra -----");

    // Create a filled vector of specified size
    let mut a = Vec64::new( 4, 0.0 ); 
    a[ 0 ] = 1.0; a[ 1 ] = 3.0; a[ 2 ] = 1.0; a[ 3 ] = - 1.0;
    println!("  * a = {}", a.clone() );
    
    // Create a vector from a standard vec
    let mut b = Vector::<f64>::create( vec![ 2.0, 1.0, 3.0, 4.0 ] );
    println!("  * b = {}", b.clone() );

    // Create an empty vector of unspecified size
    let mut c = Vec64::empty();
    println!("  * c = {}", c );

    // Do some basic algebra
    c = a.clone() + b.clone();
    println!( "  * c = a + b = {}", c );
    c = b.clone() - a.clone();
    println!( "  * c = b - a = {}", c );
    c = 3.0 * a.clone();
    println!( "  * c = 3 * a = {}", c );
    c = a.clone() / 3.0;
    println!( "  * c = a / 3 = {}", c );

    // Add some new elements and swap them around
    a.push( 2.0 );
    b.resize( 5 );
    b[ 4 ] = - 3.0;
    b.swap( 3, 4 );
    println!("  * a = {}", a );
    println!("  * b = {}", b );

    // Remove elements 
    a.pop();
    a.pop();
    b.resize( 3 );
    println!( "  * a = {}", a );
    println!( "  * b = {}", b );

    // Magnitude of the vectors 
    println!( "  * |a| = {}", a.norm_2() );
    println!( "  * |b| = {}", b.norm_p( 2.0 ) );

    // Angle between the two vectors
    let cos: f64 = a.dot( &b ) / ( a.norm_2() * b.norm_2() );
    let rad = f64::acos( cos );
    let deg = rad * 180.0 / PI;
    println!( "  * angle between a and b = {} degrees", deg );

    // Create vectors of spaced elements
    let a = Vec64::linspace( 0.0, 1.0, 5 );
    let b = Vec64::powspace( 0.0, 1.0, 5, 1.5 );
    println!( "  * Linearly spaced elements = {}", a );
    println!( "  * Power law spaced elements (p=1.5) = {}", b );
    println!( "--- FINISHED ---" );
}

