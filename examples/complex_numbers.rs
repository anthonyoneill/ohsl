// cargo run --example complex_numbers

extern crate ohsl;

pub use ohsl::complex::{Complex, Cmplx};
pub use ohsl::constant::I;

fn main() {
    println!("----- Complex numbers -----");

    // Create a new complex number
    let z = Cmplx::new( 1.0, -1.0 );
    println!( "  * z = {}", z );
    println!( "  * Re[z] = {}", z.real );
    println!( "  * Im[z] = {}", z.imag );

    // Take the conjugate
    let zbar = z.conj();
    println!("  * zbar = {}", zbar );

    // Modulus and arguement
    println!( "  * |z| = {}", z.abs() );
    println!( "  * arg(z) = {}", z.arg() );

    // Arithmetic
    let z1 = Complex::<f64>::new( 2.0, 1.0 );
    let z2 = z;
    println!( "  * z1 = {}", z1 );
    println!( "  * z2 = {}", z2 );
    println!( "  * z1 + z2 = {}", z1 + z2 );
    println!( "  * z1 - z2 = {}", z1 - z2 );
    println!( "  * z1 * z2 = {}", z1 * z2 );
    println!( "  * z1 / z2 = {}", z1 / z2 );
    let real: f64 = 2.0;
    println!( "  * z1 + 2 = {}", z1 + real );
    println!( "  * z1 - 2 = {}", z1 - real );
    println!( "  * z1 * 2 = {}", z1 * real );
    println!( "  * z1 / 2 = {}", z1 / real );
    println!( "  * z1 + i = {}", z1 + I );  
    println!( "--- FINISHED ---" );
}
