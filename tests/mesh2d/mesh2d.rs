use ohsl::{ mesh2d::Mesh2D, vector::Vec64};

#[test]
fn construction() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh[(0,0)][1] = 3.14;
    assert_eq!( mesh[(0,0)][1], 3.14 );
    assert_eq!( mesh.nvars(), 3 );
    assert_eq!( mesh.nnodes().0, 11 );
    assert_eq!( mesh.nnodes().1, 21 );
}

#[test]
fn coord_nodes() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    assert_eq!( mesh.coord( 2, 4 ), ( 0.2, 0.4 ) );
    assert!( ( mesh.xnodes()[3] - 0.3 ).abs() < 1.0e-8 );
    assert!( ( mesh.ynodes()[19] - 1.9 ).abs() < 1.0e-8 );
}

#[test]
fn set_nodes_vars() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    let vec = Vec64::new( 3, 3.14 );
    mesh.set_nodes_vars( 2, 4, vec );
    assert_eq!( mesh[(2,4)][0], 3.14 );
    assert_eq!( mesh[(2,4)][1], 3.14 );
    assert_eq!( mesh[(2,4)][2], 3.14 );
}

#[test]
fn get_nodes_vars() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh[(0,0)][1] = 3.14;
    mesh[(0,0)][2] = 6.28;
    let vec = mesh.get_nodes_vars( 0, 0 );
    assert_eq!( vec[0], 0.0 );
    assert_eq!( vec[1], 3.14 );
    assert_eq!( vec[2], 6.28 );
}

#[test]
fn assign() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh.assign( 2.718 );
    for i in 0..11 {
        for j in 0..21 {
            let vec = mesh.get_nodes_vars( i, j );
            assert_eq!( vec[0], 2.718 );
            assert_eq!( vec[1], 2.718 );
            assert_eq!( vec[2], 2.718 );
        }
    }
}

#[test]
fn cross_section() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    let vec = Vec64::new( 3, 3.14 );
    mesh.set_nodes_vars( 2, 4, vec );
    let y_section = mesh.cross_section_xnode( 2 );
    assert_eq!( y_section.nnodes(), 21 );
    assert_eq!( y_section[4][0], 3.14 );
    assert_eq!( y_section[4][1], 3.14 );
    assert_eq!( y_section[4][2], 3.14 );
    let x_section = mesh.cross_section_ynode( 4 );
    assert_eq!( x_section.nnodes(), 11 );
    assert_eq!( x_section[2][0], 3.14 );
    assert_eq!( x_section[2][1], 3.14 );
    assert_eq!( x_section[2][2], 3.14 );
}

#[test]
fn var_as_matrix() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 6 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 6 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh.assign( 2.718 );
    mesh[(0,0)][1] = 3.14;
    let mat = mesh.var_as_matrix( 1 );
    assert_eq!( mat[(0,0)], 3.14 );
    assert_eq!( mat[(0,1)], 2.718 );
}

/*#[test]
fn test_mesh2d_output() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    for i in 0..mesh.xnodes().size() {
        let x = mesh.xnodes()[i].clone();
        for j in 0..mesh.ynodes().size() {
        let y = mesh.ynodes()[j].clone();
            mesh[(i,j)][0] = 2.0 * x * y;
            mesh[(i,j)][1] = x * x + y * y;
            mesh[(i,j)][2] = x - y;
        }
    }
    //mesh.output( "./output.txt", 5 );
    mesh.output_var( "./output.txt", 0, 5)
}*/

fn two_x_y( x: f64, y: f64 ) -> f64 {
    2.0 * x * y
}

#[test]
fn apply() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh.apply( &two_x_y, 0 );
    assert!( ( mesh[(1,1)][0] - 0.02 ).abs() < 1.0e-8 );
}

#[test]
fn trapezium() {
    let x_nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let y_nodes = Vec64::linspace( 0.0, 2.0, 21 );
    let mut mesh = Mesh2D::<f64>::new( x_nodes, y_nodes, 3 );
    mesh.apply( &two_x_y, 0 );
    let integral = mesh.trapezium( 0 );
    assert!( ( integral - 2.0 ).abs() < 0.001 );
}