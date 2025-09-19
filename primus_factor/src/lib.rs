#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod shoup_factor;

pub trait LazyFactorMul<T> {
    /// Calculates `self * b (mod modulus)`.
    fn lazy_factor_mul_modulo(self, b: T, modulus: T) -> T;
}

pub trait FactorMul<T>: LazyFactorMul<T> {
    /// Calculates `self * b (mod modulus)`.
    fn factor_mul_modulo(self, b: T, modulus: T) -> T;
}

pub use shoup_factor::ShoupFactor;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use shoup_factor::SimdShoupFactor;
