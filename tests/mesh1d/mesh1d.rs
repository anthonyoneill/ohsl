use ohsl::{ mesh1d::Mesh1D, vector::Vec64};

#[test]
fn constructor_indexing() {
    let nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes, 3 );
    mesh[0][1] = 3.14;
    assert_eq!( mesh.nnodes(), 11 );
    assert_eq!( mesh.nvars(), 3 );
    assert_eq!( mesh[0][1], 3.14 );
}

#[test]
fn coord() {
    let nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let mesh = Mesh1D::<f64, f64>::new( nodes, 3 );
    assert_eq!( mesh.coord( 1 ), 0.1 );
}

#[test]
fn set_nodes_vars() {
    let nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes, 3 );
    let vec = Vec64::new( 3, 3.14 );
    mesh.set_nodes_vars( 4, vec );
    assert_eq!( mesh[4][0], 3.14 );
    assert_eq!( mesh[4][1], 3.14 );
    assert_eq!( mesh[4][2], 3.14 );
}

#[test]
fn get_nodes_vars() {
    let nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes, 3 );
    let vec = Vec64::new( 3, 3.14 );
    mesh.set_nodes_vars( 4, vec );
    assert_eq!( mesh.nodes()[10], 1.0 );
}

#[test]
fn get_interpolated_vars() {
    let nodes = Vec64::linspace( 0.0, 1.0, 11 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes.clone(), 1 );
    for i in 0..nodes.size() {
        mesh[i][0] = 2.0 * nodes[i].clone();
    }
    let vec = mesh.get_interpolated_vars( 0.55 );
    assert_eq!( vec[0], 1.1 );
}

/*#[test]
fn output() {
    let mut nodes = Vec64::empty();
    nodes.linspace( 0.0, 1.0, 11 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes.clone(), 3 );
    for i in 0..nodes.size() {
        let x = nodes[i].clone();
        mesh[i][0] = 2.0 * x;
        mesh[i][1] = x * x;
        mesh[i][2] = 0.5 * x;
    }
    mesh.output( "./output.txt", 5 );
}*/

#[test]
fn trapezium() {
    let nodes = Vec64::linspace( 0.0, 1.0, 21 );
    let mut mesh = Mesh1D::<f64, f64>::new( nodes.clone(), 2 );
    for i in 0..nodes.size() {
        let x = nodes[i].clone();
        mesh[i][0] = 2.0 * x;
        mesh[i][1] = x * x;
    }
    assert!( ( mesh.trapezium( 0 ) - 1.0 ).abs() < 0.001 );
    assert!( ( mesh.trapezium( 1 ) - 1.0 / 3.0).abs() < 0.001 );
}

/*#[test]
fn read() {
    let nodes = Vec64::empty();
    let mut mesh = Mesh1D::<f64, f64>::new( nodes, 3 );
    mesh.read( "./output.txt" );
    assert_eq!( mesh[2][0], 0.4 );
    assert_eq!( mesh[8][1], 0.64 );
    assert_eq!( mesh[10][2], 0.5 );
    assert_eq!( mesh.nnodes(), 11 );
}*/