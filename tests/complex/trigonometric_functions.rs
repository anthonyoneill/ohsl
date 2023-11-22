use ohsl::complex::Cmplx;

#[test]
fn test_complex_sin() {
    let z = Cmplx::new( 1.0, 1.0 );
    let s = z.sin();
    assert!( (s - Cmplx::new( 1.2984575814159773, 0.6349639147847361 )).abs() < 1e-15 );
}

#[test]
fn test_complex_cos() {
    let z = Cmplx::new( 1.0, 1.0 );
    let c = z.cos();
    assert!( (c - Cmplx::new( 0.8337300251311491, -0.9888977057628651 )).abs() < 1e-15 );
}

#[test]
fn test_complex_tan() {
    let z = Cmplx::new( 1.0, 1.0 );
    let t = z.tan();
    assert!( (t - Cmplx::new( 0.2717525853195117, 1.0839233273386946 )).abs() < 1e-15 );
}

#[test]
fn test_complex_sec() {
    let z = Cmplx::new( 1.0, 1.0 );
    let s = z.sec();
    assert!( (s - Cmplx::new( 0.4983370305551868, 0.591083841721045 )).abs() < 1e-15 );
}

#[test]
fn test_complex_csc() {
    let z = Cmplx::new( 1.0, 1.0 );
    let c = z.csc();
    assert!( (c - Cmplx::new( 0.621518017170428, -0.30393100162842697 )).abs() < 1e-15 );
}

#[test]
fn test_complex_cot() {
    let z = Cmplx::new( 1.0, 1.0 );
    let c = z.cot();
    assert!( (c - Cmplx::new( 0.21762156185440243, -0.8680141428959249 )).abs() < 1e-15 );
}

#[test]
fn test_complex_asin() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.asin();
    assert!( (a - Cmplx::new( 0.6662394324925153, 1.0612750619050355 )).abs() < 1e-15 );
}

#[test]
fn test_complex_acos() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.acos();
    assert!( (a - Cmplx::new( 0.9045568943023813, -1.0612750619050355 )).abs() < 1e-15 );
}

#[test]
fn test_complex_atan() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.atan();
    assert!( (a - Cmplx::new( 1.0172219678978514, 0.4023594781085251 )).abs() < 1e-15 );
}

#[test]
fn test_complex_asec() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.asec();
    assert!( (a - Cmplx::new( 1.1185178796437059, 0.5306375309525178 )).abs() < 1e-15 );
}

#[test]
fn test_complex_acsc() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.acsc();
    assert!( (a - Cmplx::new( 0.45227844715119065, -0.5306375309525178 )).abs() < 1e-15 );
}

#[test]
fn test_complex_acot() {
    let z = Cmplx::new( 1.0, 1.0 );
    let a = z.acot();
    assert!( (a - Cmplx::new( 0.5535743588970451, -0.4023594781085251 )).abs() < 1e-15 );
}
