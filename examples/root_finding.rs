// cargo run --example root_finding

extern crate ohsl;

pub use ohsl::{vector::Vec64, matrix::Mat64, newton::Newton};

fn main() {
    println!( "------------------ Root finding ------------------" );

    println!( " * Let us solve the equation cos(x) - x = 0 ");
    println!( " * with initial guess x_0 = 1.");

    fn function(x: f64) -> f64 {
        libm::cos( x ) - x
    }
    let newton = Newton::<f64>::new( 1.0 );
    let solution = newton.solve( &function );
    println!( " * Our solution is x = {:.6}", solution.unwrap() );
    println!( " * to six decimal places.");
    println!( "---------------------------------------------------" );
    println!( " * Now we shall solve the set of equations ");
    println!( " * x^3 + y - 1 = 0, " );
    println!( " * y^3 - x + 1 = 0, " );
    println!( " * using Newton's method and the initial guess " );
    let guess = Vec64::create( vec![ 0.5, 0.25 ] );
    println!( " * ( x_0, y_0 ) = {}.", guess.clone() );

    fn vector_function( x: Vec64 ) -> Vec64 {
        let mut f = Vec64::new( 2, 0.0 );
        f[0] = libm::pow( x[0], 3.0 ) + x[1] -1.0;
        f[1] = libm::pow( x[1], 3.0 ) - x[0] + 1.0;
        f
    }
    let newton = Newton::<Vec64>::new( guess );
    let solution = newton.solve( &vector_function ).unwrap();
    println!( " * Our solution is x = {:.2} and", solution[0] );
    println!( " * y = {:.2} to two decimal places.", solution[1] );
    println!( "-------------------- FINISHED ---------------------" );
}

