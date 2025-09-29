#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod integer_traits;

mod integer;
mod unsigned_integer;

mod big_integer;

#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

pub use integer_traits::*;
pub use primus_utils::ByteCount;

pub use integer::Integer;
pub use unsigned_integer::UnsignedInteger;

pub use big_integer::*;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::{SimdArray, SimdMaskArray, SimdUnsignedInteger};
