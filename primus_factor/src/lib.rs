#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod shoup_factor;

pub trait LazyFactorMul<T, M> {
    /// Calculates `self * b (mod modulus)`.
    fn lazy_factor_mul_modulo(self, b: T, modulus: M) -> T;
}

pub trait FactorMul<T, M> {
    /// Calculates `self * b (mod modulus)`.
    fn factor_mul_modulo(self, b: T, modulus: M) -> T;
}

pub use shoup_factor::ShoupFactor;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use shoup_factor::SimdShoupFactor;
