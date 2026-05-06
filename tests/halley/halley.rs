use ohsl::halley::Halley;
use ohsl::complex::Cmplx;
use ohsl::vector::Vec64;

#[test]
fn constructor() {
    let halley = Halley::<f64>::new( 0.0 );
    let parameters = halley.parameters();
    assert_eq!( parameters.0, 1.0e-8 );
    assert_eq!( parameters.1, 1.0e-6 );
    assert_eq!( parameters.2, 20 );
    assert_eq!( parameters.3, 0.0 );
}

#[test]
fn parameters() {
    let mut halley = Halley::<f64>::new( 0.0 );
    halley.tol = 1.0e-6;
    halley.delta = 1.0e-7;
    halley.max_iter = 25;
    halley.guess = 1.0;
    let parameters = halley.parameters();
    assert_eq!( parameters.0, 1.0e-6 );
    assert_eq!( parameters.1, 1.0e-7 );
    assert_eq!( parameters.2, 25 );
    assert_eq!( parameters.3, 1.0 );
}

fn x_squared_minus_four(x: &f64) -> f64 {
    x * x - 4.0
}

#[test]
fn solve_f64() {
    let mut halley = Halley::<f64>::new( 1.0 );
    let solution = halley.solve( &x_squared_minus_four );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, 2.0 );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = -1.0;
    let solution = halley.solve( &x_squared_minus_four );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, -2.0 );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
}

fn x_squared_minus_four_first(x: &f64) -> f64 {
    2.0 * x
}

fn x_squared_minus_four_second(_x: &f64) -> f64 {
    2.0
}

#[test]
fn solve_derivative_f64() {
    let mut halley = Halley::<f64>::new( 1.0 );
    let solution = halley.solve_derivative( 
        &x_squared_minus_four, 
        &x_squared_minus_four_first, 
        &x_squared_minus_four_second 
    );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, 2.0 );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = -1.0;
    let solution = halley.solve_derivative( 
        &x_squared_minus_four, 
        &x_squared_minus_four_first, 
        &x_squared_minus_four_second 
    );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, -2.0 );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
}

fn x_cubed_minus_two(x: &Cmplx) -> Cmplx {
    x * x * x - 2.0
}

#[test]
fn solve_cmplx() {
    let mut halley = Halley::<Cmplx>::new( Cmplx::new( 1.0, 0.0 ) );
    let solution = halley.solve( &x_cubed_minus_two );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( 1.2599210498948732, 0.0 ) );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = Cmplx::new( -1.0, 1.0 );
    let solution = halley.solve( &x_cubed_minus_two );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( -0.6299605249474366, 1.0911236359717214 ) );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = Cmplx::new( -1.0, -1.0 );
    let solution = halley.solve( &x_cubed_minus_two );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( -0.6299605249474366, -1.0911236359717214 ) );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
}

fn x_cubed_minus_two_first(x: &Cmplx) -> Cmplx {
    3.0 * x * x
}

fn x_cubed_minus_two_second(x: &Cmplx) -> Cmplx {
    6.0 * x
}

#[test]
fn solve_cmplx_derivative() {
    let mut halley = Halley::<Cmplx>::new( Cmplx::new( 1.0, 0.0 ) );
    let solution = halley.solve_derivative( 
        &x_cubed_minus_two, &x_cubed_minus_two_first, &x_cubed_minus_two_second 
    );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( 1.2599210498948732, 0.0 ) );
            println!( "Halley with derivative converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley with derivative failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = Cmplx::new( -1.0, 1.0 );
    let solution = halley.solve_derivative( 
        &x_cubed_minus_two, &x_cubed_minus_two_first, &x_cubed_minus_two_second 
    );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( -0.6299605249474366, 1.0911236359717214 ) );
            println!( "Halley with derivative converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley with derivative failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
    halley.guess = Cmplx::new( -1.0, -1.0 );
    let solution = halley.solve_derivative( 
        &x_cubed_minus_two, &x_cubed_minus_two_first, &x_cubed_minus_two_second 
    );
    match solution {
        Ok( (sol, iter) ) => {
            assert_eq!( sol, Cmplx::new( -0.6299605249474366, -1.0911236359717214 ) );
            println!( "Halley with derivative converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley with derivative failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
}

fn vecfunc( x: &Vec64 ) -> Vec64 {
    let mut f = Vec64::new( 2, 0.0 );
    f[0] = f64::powf( x[0], 3.0 ) + x[1] -1.0;
    f[1] = f64::powf( x[1], 3.0 ) - x[0] + 1.0;
    /*
        x^3 + y - 1 = 0,
        y^3 - x + 1 = 0,
        (x,y) = (1,0) is the only (real) solution
    */
    f
}

#[test]
fn solve_vec64() {
    let guess = Vec64::create( vec![ 0.5, 0.25 ] );
    let halley = Halley::<Vec64>::new( guess );
    let solution = halley.solve( &vecfunc );
    match solution {
        Ok( (sol, iter) ) => {
            assert!( ( sol[0] - 1.0).abs() < 1.0e-8 );
            assert!( ( sol[1] - 0.0).abs() < 1.0e-8 );
            println!( "Halley converged in {} iterations", iter );
        },
        Err( sol ) => {
            println!( "Halley failed to converge, last solution was {}", sol );
            assert!( false );
        }
    }
}