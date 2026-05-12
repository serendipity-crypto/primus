#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

mod multiply;
mod shoup_factor;

pub trait LazyFactorMul<T> {
    /// Calculates `self * b (mod modulus)`.
    fn lazy_factor_mul_modulo(self, b: T, modulus: T) -> T;
}

pub trait FactorMul<T>: LazyFactorMul<T> {
    /// Calculates `self * b (mod modulus)`.
    fn factor_mul_modulo(self, b: T, modulus: T) -> T;
}

/// Slice-level multiplication by a precomputed factor.
///
/// Implementations may use SIMD internally when the `nightly` and `simd`
/// features are enabled. Callers keep the normal scalar slice layout and the
/// remainder is handled by the scalar path.
pub trait LazyFactorSliceOps<T> {
    /// Calculates `factor * value (mod 2*modulus)` for each element in-place.
    fn lazy_factor_mul_slice_assign(self, values: &mut [T], modulus: T);

    /// Calculates `factor * input (mod 2*modulus)` into `output`.
    ///
    /// # Panics
    ///
    /// Panics if `input.len() != output.len()`.
    fn lazy_factor_mul_slice_to(self, input: &[T], output: &mut [T], modulus: T);
}

/// Slice-level multiplication by a precomputed factor.
///
/// Implementations may use SIMD internally when the `nightly` and `simd`
/// features are enabled. Callers keep the normal scalar slice layout and the
/// remainder is handled by the scalar path.
pub trait FactorSliceOps<T> {
    /// Calculates `factor * value (mod modulus)` for each element in-place.
    fn factor_mul_slice_assign(self, values: &mut [T], modulus: T);

    /// Calculates `factor * input (mod modulus)` into `output`.
    ///
    /// # Panics
    ///
    /// Panics if `input.len() != output.len()`.
    fn factor_mul_slice_to(self, input: &[T], output: &mut [T], modulus: T);

    /// Calculates `acc += factor * rhs (mod modulus)` element-wise.
    ///
    /// # Panics
    ///
    /// Panics if `acc.len() != rhs.len()`.
    fn add_factor_mul_slice_assign(self, acc: &mut [T], rhs: &[T], modulus: T);
}

pub use multiply::MultiplyFactor;
pub use shoup_factor::ShoupFactor;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use shoup_factor::{SimdShoupFactor, default_lanes, simd_kernel};
