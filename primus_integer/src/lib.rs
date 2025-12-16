#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod macros;

mod integer_traits;

mod integer;
mod unsigned_integer;

mod big_integer;

mod data;

#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

pub mod size;

pub use integer_traits::*;

pub use integer::Integer;
pub use unsigned_integer::UnsignedInteger;

pub use big_integer::*;

pub use data::{Data, DataMut, DataOwned, RawData};

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::{SimdArray, SimdMaskArray, SimdUnsignedInteger};
