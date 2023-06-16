use ohsl::vector::Vector;
use ohsl::complex::Cmplx;

#[test]
fn conjugate() {
    let v = Vector::<Cmplx>::new( 5, Cmplx::new( 1.0, 2.0 ) );
    assert_eq!( v[0].real, 1.0 );
    let conj = v.conj();
    assert_eq!( conj[0].imag, -2.0 );
}

#[test]
fn real() {
    let v = Vector::<Cmplx>::new( 5, Cmplx::new( 1.0, 2.0 ) );
    let real = v.real();
    assert_eq!( real[0], 1.0 );
}