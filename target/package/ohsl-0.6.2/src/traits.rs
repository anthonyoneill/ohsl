use core::ops::{Add, Div, Mul, Neg, Sub};
use core::ops::{AddAssign, DivAssign, MulAssign, SubAssign};

pub trait Operations<Rhs = Self, Output = Self>:
      Add<Rhs, Output = Output>
    + Sub<Rhs, Output = Output>
    + Mul<Rhs, Output = Output>
    + Div<Rhs, Output = Output>
{
}

impl<T, Rhs, Output> Operations<Rhs, Output> for T where
    T:    Add<Rhs, Output = Output>
        + Sub<Rhs, Output = Output>
        + Mul<Rhs, Output = Output>
        + Div<Rhs, Output = Output>
{
}

pub trait AssignOperations<Rhs = Self>:
      AddAssign<Rhs> 
    + SubAssign<Rhs> 
    + MulAssign<Rhs> 
    + DivAssign<Rhs> 
{
}

impl<T, Rhs> AssignOperations<Rhs> for T where 
    T:    AddAssign<Rhs> 
        + SubAssign<Rhs> 
        + MulAssign<Rhs> 
        + DivAssign<Rhs>
{
}

pub trait Zero: Sized + Add<Self, Output = Self> {
    /// Returns the additive identity element 0
    fn zero() -> Self;
}

pub trait One: Sized + Mul<Self, Output = Self> {
    /// Returns the multiplicative identity element 1
    fn one() -> Self;
}
   

pub trait Number: PartialEq + Sized + Operations + AssignOperations + Zero + One
{

}

pub trait Signed: Number + Neg<Output = Self> {
    fn abs(&self) -> Self;
}

macro_rules! impl_identity {
    ($name: ident for $($t: ty)*, $method: ident, $v: expr) => ($(
        impl $name for $t {
            #[inline]
            fn $method() -> $t {
                $v
            }
        }
    )*)
}

impl_identity!( Zero for f32 f64, zero, 0.0 );
impl_identity!( One for f32 f64, one, 1.0 );
impl_identity!( Zero for usize u8 u16 u32 u64 isize i8 i16 i32 i64, zero, 0 );
impl_identity!( One for usize u8 u16 u32 u64 isize i8 i16 i32 i64, one, 1 );

macro_rules! impl_signed {
    ($($t:ty)*, $v: expr) => ($(
        impl Signed for $t {
            #[inline]
            fn abs(&self) -> $t {
                if *self < $v { -*self } else { *self }
            }
        }
    )*)
}

impl_signed!( isize i8 i16 i32 i64, 0 );
impl_signed!( f32 f64, 0.0 );

macro_rules! impl_trait {
    ($name: ident for $($t: ty)*) => ($(
        impl $name for $t {

        }
    )*)
}

impl_trait!( Number for usize u8 u16 u32 u64 isize i8 i16 i32 i64 f32 f64);