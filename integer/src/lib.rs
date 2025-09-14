#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod integer_traits;

mod integer;
mod unsigned_integer;

#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

pub use integer_traits::*;

pub use integer::Integer;
pub use unsigned_integer::UnsignedInteger;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::{SimdArray, SimdMaskArray, SimdUnsignedInteger};
