// cargo run --example linear_systems

extern crate ohsl;

pub use ohsl::vector::{Vector, Vec64};
pub use ohsl::matrix::{Matrix, Mat64};
pub use ohsl::complex::{Complex, Cmplx};



fn main() {
    println!( "------------------ Linear system ------------------" );
    println!( "  * Solving the linear system Ax = b, for x, where" );
    let mut a = Mat64::new( 3, 3, 0.0 );
    a[(0,0)] = 1.0; a[(0,1)] = 1.0; a[(0,2)] =   1.0;
    a[(1,0)] = 0.0; a[(1,1)] = 2.0; a[(1,2)] =   5.0;
    a[(2,0)] = 2.0; a[(2,1)] = 5.0; a[(2,2)] = - 1.0;
    println!( "  * A ={}", a );
    let b = Vec64::create( vec![ 6.0, -4.0, 27.0 ] ); 
    println!( "  * b^T ={}", b );
    let x = a.clone().solve_basic( &b );
    println!( "  * gives the solution vector");
    println!( "  * x^T ={}", x );
    println!( "---------------------------------------------------" );
    println!( "  * Lets solve the complex linear system Cy = d where" );
    let mut c = Matrix::<Cmplx>::new( 2, 2, Cmplx::new( 1.0, 1.0 ) );
    c[(0,1)] = Cmplx::new( -1.0, 0.0 );
    c[(1,0)] = Cmplx::new( 1.0, -1.0 );
    println!( "  * C ={}", c );
    let mut d = Vector::<Cmplx>::new( 2, Cmplx::new( 0.0, 1.0 ) );
    d[1] = Cmplx::new( 1.0, 0.0 );
    println!( "  * d^T ={}", d );
    println!( "  * Using LU decomposition we find that ");
    let y = c.clone().solve_lu( &d );
    println!( "  * y^T ={}", y );
    println!( "---------------------------------------------------" );
    println!( "  * We may also find the determinant of a matrix" );
    println!( "  * |A| = {}", a.determinant() );
    println!( "  * |C| = {}", c.determinant() ); 
    println!( "  * or the inverse of a matrix" );
    let inverse = a.inverse();
    println!( "  * A^-1 =\n{}", inverse );
    println!( "  * C^-1 =\n{}", c.inverse() );
    println!( "  * We can check this by multiplication" );
    println!( "  * I = A * A^-1 =\n{}", a * inverse );
    println!( "  * which is sufficiently accurate for our purposes.");
    println!( "-------------------- FINISHED ---------------------" );
}

