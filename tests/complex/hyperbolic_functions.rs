use ohsl::complex::Complex;

#[test]
fn test_complex_sinh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.sinh() - Complex::new( 0.6349639147847361, 1.2984575814159773 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_cosh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.cosh() - Complex::new( 0.8337300251311491, 0.9888977057628651 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_tanh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.tanh() - Complex::new( 1.0839233273386944, 0.27175258531951174 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_sech() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.sech() - Complex::new( 0.49833703055518686, -0.5910838417210451 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_csch() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.csch() - Complex::new( 0.3039310016284264, -0.6215180171704284 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_coth() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.coth() - Complex::new( 0.8680141428959249, -0.21762156185440212 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_asinh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.asinh() - Complex::new( 1.0612750619050355, 0.6662394324925153 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_acosh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.acosh() - Complex::new( 1.0612750619050355, 0.9045568943023813 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_atanh() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.atanh() - Complex::new( 0.4023594781085251, 1.0172219678978514 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_asech() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.asech() - Complex::new( 0.5306375309525178, -1.118517879643705 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_acsch() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.acsch() - Complex::new( 0.5306375309525178, -0.45227844715119065 ) ).abs() < 1e-15);
}

#[test]
fn test_complex_acoth() {
    let z = Complex::new(1.0, 1.0);
    assert!((z.acoth() - Complex::new( 0.4023594781085251, -0.5535743588970452 ) ).abs() < 1e-15);
}
