#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod shoup_factor;

pub trait LazyMulFactor<T, M> {
    /// Calculates `self * b (mod modulus)`.
    fn lazy_mul_modulo(self, b: T, modulus: M) -> T;
}

pub trait MulFactor<T, M> {
    /// Calculates `self * b (mod modulus)`.
    fn mul_modulo(self, b: T, modulus: M) -> T;
}

pub use shoup_factor::ShoupFactor;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use shoup_factor::SimdShoupFactor;
