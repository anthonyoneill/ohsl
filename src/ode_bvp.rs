pub use crate::vector::{Vector, Vec64};
pub use crate::mesh1d::Mesh1D;
pub use crate::banded::Banded;
pub use crate::matrix::Mat64;

pub use crate::traits::Number;

pub struct ODEBVP<T> {
    pub max_iter: usize,            // Maximum number of iterations
    pub tol: f64,                   // Tolerance for convergence
    pub delta: f64,                 // Finite difference step size
    pub solution: Mesh1D<T, f64>,   // Solution mesh
}   

impl<T: Clone + Number> ODEBVP<T> {
    /// Create a new ODEBVP solver 
    #[inline]
    pub fn new( nodes: Vector<f64>, order: usize ) -> Self {
        let max_iter: usize = 20;
        let tol: f64 = 1.0e-8;
        let delta: f64 = 1.0e-8;
        let solution = Mesh1D::<T, f64>::new( nodes, order );
        ODEBVP { max_iter, tol, delta, solution }
    }

    /// Edit the maximum number of iterations 
    #[inline]
    pub fn iterations(&mut self, iterations: usize ) {
        self.max_iter = iterations;
    }

    /// Edit the convergence tolerance 
    #[inline]
    pub fn tolerance(&mut self, tolerance: f64 ) {
        self.tol = tolerance;
    }

    /// Edit the finite difference derivative step 
    #[inline]
    pub fn delta(&mut self, delta: f64 ) {
        self.delta = delta;
    }
}

impl ODEBVP<f64> {
    pub fn solve(&mut self, 
        eqn: &dyn Fn(&Vec64, &f64) -> Vec64, 
        bc_left: &dyn Fn(&Vec64) -> Vector<Option<f64>>,
        bc_right: &dyn Fn(&Vec64) -> Vector<Option<f64>>,
        jacobian: Option<&dyn Fn(&Vec64, &f64) -> Mat64>
    ) -> Result<(), f64> {
        //TODO should probably check sizes of eqn, bc_left, bc_right here
        let order = self.solution.nvars();
        let n = self.solution.nodes().size();
        let mut max_residual = 1.0;

        let size = n * order;
        let off_diag = 2 * order - 1;
        let mut a_mat = Banded::<f64>::new( size, off_diag, off_diag, 0.0 );
        let mut b_vec = Vec64::new( size, 0.0 );

        for _iter in 0..self.max_iter {
            self.assemble_matrix_problem( eqn, bc_left, bc_right, jacobian, &mut a_mat, &mut b_vec );
            let sol = a_mat.solve( &b_vec );
            // Update the solution mesh
            for i in 0..n {
                for j in 0..order {
                    self.solution[i][j] += sol[ i * order + j ];
                }
                
            }
            max_residual = b_vec.norm_inf();
            //println!( "max residual: {:.2e}", max_residual );
            
            if max_residual <= self.tol {
                return Ok( () )
            }
        }


        Err( max_residual )
    }

    fn assemble_matrix_problem( &mut self, 
        eqn: &dyn Fn(&Vec64, &f64) -> Vec64, 
        bc_left: &dyn Fn(&Vec64) -> Vector<Option<f64>>, 
        bc_right: &dyn Fn(&Vec64) -> Vector<Option<f64>>,
        jacobian: Option<&dyn Fn(&Vec64, &f64) -> Mat64>,
        a: &mut Banded<f64>,
        b: &mut Vec64 
    ) {
        let order = self.solution.nvars();
        let n = self.solution.nodes().size();

        a.fill( 0.0 );          // Clear the matrix
        b.assign( 0.0 );        // Clear the rhs vector
        let mut row = 0;        // Row counter

        // Left boundary condition
        let bc_left_guess: Vec64 = self.solution.get_nodes_vars( 0 );
        let bc_left_residual: Vector<Option<f64>> = bc_left( &bc_left_guess );

        for i in 0..order {
            if let Some( res ) = bc_left_residual[ i ] {
                a[ ( row, i ) ] = 1.0;
                b[ row ] = -res;
                row += 1;
            }
        }

        // Interior nodes
        for i in 0..n-1 {
            let x_i = self.solution.coord( i );
            let x_i_1 = self.solution.coord( i + 1 );
            let dx_i = x_i_1 - x_i;
            let inv_dx_i = 1.0 / dx_i;
            let x_i_half = 0.5 * ( x_i + x_i_1 );

            let y_g_i = self.solution.get_nodes_vars( i );
            let y_g_i_1 = self.solution.get_nodes_vars( i + 1 );
            let y_g_i_half = 0.5 * ( &y_g_i + &y_g_i_1 );

            // Function evaluation at the midpoint i + 1/2
            let f_fun = eqn( &y_g_i_half, &x_i_half );

            let mut jac = Mat64::new( order, order, 0.0 );

            if let Some(jacobian_func) = jacobian {
                // Use provided Jacobian function
                jac = jacobian_func( &y_g_i_half, &x_i_half );
            } else {
                // Approximate Jacobian at the midpoint i + 1/2 via finite differences
                for alpha in 0..order {
                    for beta in 0..order {
                        let mut delta = Vec64::new( order, 0.0 );
                        delta[ beta ] = self.delta;
                        let y_g_i_half_star = &y_g_i_half + &delta;

                        // Perturbed function evaluation at i + 1/2
                        let f_fun_star = eqn( &y_g_i_half_star, &x_i_half );

                        // Approximate Jacobian entry
                        jac[( alpha, beta )] = ( f_fun_star[ alpha ] - f_fun[ alpha ] ) / self.delta;
                    }
                }
            }

            let dy_g_i = inv_dx_i * ( &y_g_i_1 - &y_g_i );
            let r = &f_fun - &dy_g_i;

            // Coefficient matrices
            let mut identity = Mat64::new( order, order, 0.0 );
            identity.fill_diag( inv_dx_i );
            let jac_half = &jac * 0.5;
            let k = &identity - &jac_half;
            let l = &-identity - &jac_half;

            // Fill the rows in the matrix and rhs vector
            for n in 0..order {
                for j in 0..order {
                    // Left node contribution
                    a[ ( row, ( i * order ) + j ) ] = l[( n, j )];
                    // Right node contribution
                    a[ ( row, ( ( i + 1 ) * order ) + j ) ] = k[( n, j )];
                }
                b[ row ] = r[ n ];
                row += 1;
            }
            
        }

        // Right boundary conditions
        let bc_right_guess: Vec64 = self.solution.get_nodes_vars( n - 1 );
        let bc_right_residual: Vector<Option<f64>> = bc_right( &bc_right_guess );

        for i in 0..order {
            if let Some( res ) = bc_right_residual[ i ] {
                a[ ( row, ( ( n - 1 ) * order ) + i ) ] = 1.0;
                b[ row ] = -res;
                row += 1;
            }
        }
    }

}