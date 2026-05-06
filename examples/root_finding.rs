// cargo run --example root_finding

extern crate ohsl;

pub use ohsl::{
    vector::Vec64, 
    matrix::Mat64, 
    newton::Newton,
    halley::Halley
};

fn main() {
    println!( "------------------ Root finding ------------------" );
    println!( " * Let us solve the equation cos(x) - x = 0 ");
    println!( " * with initial guess x_0 = 1.");

    fn function(x: &f64) -> f64 {
        f64::cos( *x ) - *x
    }
    let newton = Newton::<f64>::new( 1.0 );
    let (solution, iterations) = newton.solve( &function ).unwrap();
    println!( " * Our solution is x = {:.6} to six decimal places.", solution );
    println!( " * Newton's method converged in {} iterations.", iterations );
    println!( "-------------------------------------------------------------" );
    println!( " * Now we shall solve the set of equations ");
    println!( " * x^3 + y - 1 = 0, " );
    println!( " * y^3 - x + 1 = 0, " );
    println!( " * using Newton's method and the initial guess " );
    let guess = Vec64::create( vec![ 0.5, 0.25 ] );
    println!( " * ( x_0, y_0 ) = {}.", guess.clone() );

    fn vector_function( x: &Vec64 ) -> Vec64 {
        let mut f = Vec64::new( 2, 0.0 );
        f[0] = f64::powf( x[0], 3.0 ) + x[1] -1.0;
        f[1] = f64::powf( x[1], 3.0 ) - x[0] + 1.0;
        f
    }
    let newton = Newton::<Vec64>::new( guess );
    let (solution, iterations) = newton.solve( &vector_function ).unwrap();
    println!( " * Our solution is x = {:.2} and", solution[0] );
    println!( " * y = {:.2} to two decimal places.", solution[1] );
    println!( " * Newton's method converged in {} iterations.", iterations );
    println!( "-------------------------------------------------------------" );
    println!( " * Let us solve the equation x^3 - x^2 - 8x + 12 = 0 ");
    println!( " * with initial guess x_0 = 1.");
    fn slow(x: &f64) -> f64 {
        let x_squared = x * x;
        x * x_squared - x_squared - 8.0 * x + 12.0
    }
    let mut newton = Newton::<f64>::new( 1.0 );
    newton.max_iter = 30;
    let (solution, iterations) = newton.solve( &slow ).unwrap();
    println!( " * Our solution is x = {:.6} to six decimal places.", solution );
    println!( " * Newton's method converged in {} iterations.", iterations );
    println!( " * Newton's method is slow for this problem." );
    println!( " * Now let's solve the same equation using Halley's method." );
    let halley = Halley::<f64>::new( 1.0 );
    let (solution, iterations) = halley.solve( &slow ).unwrap();
    println!( " * Our solution is x = {:.6} to six decimal places.", solution );
    println!( " * Halley's method converged in {} iterations.", iterations );
    println!( "-------------------------------------------------------------" );
    println!( " * Now we shall solve the set of equations ");
    println!( " * e^(-x + y) - 0.1 = 0, " );
    println!( " * e^(-x - y) - 0.1 = 0, " );
    println!( " * using Newton's method and the initial guess " );
    let guess = Vec64::create( vec![ 4.3, 2.0 ] );
    println!( " * ( x_0, y_0 ) = {}.", guess.clone() );
    fn slow_vector_function( x: &Vec64 ) -> Vec64 {
        let mut f = Vec64::new( 2, 0.0 );
        f[0] = f64::exp( -x[0] + x[1]) - 0.1;
        f[1] = f64::exp( -x[0] - x[1] ) -0.1;
        f
    }
    let mut newton = Newton::<Vec64>::new( guess.clone() );
    newton.max_iter = 100;
    let (sol, iter) = newton.solve( &slow_vector_function ).unwrap();
    println!( " * Our solution is x = {:.2} and y = {:.2} to two decimal places.", sol[0], sol[1] );
    println!( " * Newton's method converged in {} iterations.", iter );
    println!( " * Newton's method is slow for this problem." );
    println!( " * Now let's solve the same equation using Halley's method." );
    let halley = Halley::<Vec64>::new( guess );
    let (sol, iter) = halley.solve( &slow_vector_function ).unwrap();
    println!( " * Our solution is x = {:.2} and y = {:.2} to two decimal places.", sol[0], sol[1] );
    println!( " * Halley's method converged in {} iterations.", iter );
    println!( "------------------------- FINISHED --------------------------" );
}

