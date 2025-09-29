use primus_integer::{DivRemScalar, UnsignedInteger};

use crate::{FactorMul, LazyFactorMul};

#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::SimdShoupFactor;

/// A number used for fast modular multiplication.
///
/// This is efficient if many operations are multiplied by
/// the same number and then reduced with the same modulus.
#[derive(Debug, Clone, Copy, Default)]
pub struct ShoupFactor<T: UnsignedInteger> {
    /// value
    value: T,
    /// quotient
    quotient: T,
}

impl<T: UnsignedInteger> ShoupFactor<T> {
    /// Constructs a [`ShoupFactor<T>`].
    ///
    /// * `value` must be less than `modulus`.
    #[inline]
    pub fn new(value: T, modulus: T) -> Self {
        debug_assert!(value < modulus);

        // Calculate the quotient of `value * 2^64 / modulus`.
        let mut quotient = [T::ZERO; 2];
        DivRemScalar::div_rem_scalar(&[T::ZERO, value], modulus, &mut quotient);

        Self {
            value,
            quotient: quotient[0],
        }
    }

    /// Resets the `modulus` of [`ShoupFactor<T>`].
    #[inline]
    pub fn set_modulus(&mut self, modulus: T) {
        debug_assert!(self.value < modulus);

        // Calculate the quotient of `value * 2^64 / modulus`.
        let mut quotient = [T::ZERO; 2];
        DivRemScalar::div_rem_scalar(&[T::ZERO, self.value], modulus, &mut quotient);

        self.quotient = quotient[0];
    }

    /// Resets the content of [`ShoupFactor<T>`].
    ///
    /// * `value` must be less than `modulus`.
    #[inline]
    pub fn set(&mut self, value: T, modulus: T) {
        self.value = value;
        self.set_modulus(modulus);
    }

    /// Returns the value of this [`ShoupFactor<T>`].
    #[inline]
    pub const fn value(self) -> T {
        self.value
    }

    /// Returns the quotient of this [`ShoupFactor<T>`].
    #[inline]
    pub const fn quotient(self) -> T {
        self.quotient
    }
}

impl<T: UnsignedInteger> LazyFactorMul<T> for ShoupFactor<T> {
    /// Calculates `a * b mod modulus`.
    ///
    /// The result is in [0, 2 * `modulus`).
    ///
    /// # Proof
    ///
    /// Let `x = b`, `w = a.value`, `w' = a.quotient`, `p = modulus` and `β = 2^(64)`.
    ///
    /// By definition, `w' = ⌊wβ/p⌋`. Let `q = ⌊w'x/β⌋`.
    ///
    /// Then, `0 ≤ wβ/p - w' < 1`, `0 ≤ w'x/β - q < 1`.
    ///
    /// Multiplying by `xp/β` and `p` respectively, and adding, yields
    ///
    /// `0 ≤ wx - qp < xp/β + p < 2p < β`
    #[inline]
    fn lazy_factor_mul_modulo(self, b: T, modulus: T) -> T {
        let hw = self.quotient.widening_mul_hw(b);
        self.value
            .wrapping_mul(b)
            .wrapping_sub(modulus.wrapping_mul(hw))
    }
}

impl<T: UnsignedInteger> FactorMul<T> for ShoupFactor<T> {
    /// Calculates `self * b mod modulus`.
    ///
    /// The result is in [0, `modulus`).
    #[inline]
    fn factor_mul_modulo(self, b: T, modulus: T) -> T {
        let t = self.lazy_factor_mul_modulo(b, modulus);
        if t >= modulus { t - modulus } else { t }
    }
}

#[cfg(test)]
mod tests {
    use rand::{
        Rng,
        distr::{Distribution, Uniform},
    };

    use super::*;

    type ValueT = u32;
    type WideT = u64;

    const N: usize = 32;

    #[test]
    fn test_shoup() {
        let mut rng = rand::rng();

        let modulus: ValueT = 132120577;
        let distr = Uniform::new(0, modulus).unwrap();

        let shoup = ShoupFactor::new(rng.sample(distr), modulus);
        let data: Vec<ValueT> = distr.sample_iter(&mut rng).take(N).collect();

        let shoup_res: Vec<ValueT> = data
            .iter()
            .map(|&v| shoup.factor_mul_modulo(v, modulus))
            .collect();
        let normal_res: Vec<ValueT> = data
            .iter()
            .map(|&v| (v as WideT) * (shoup.value as WideT) % (modulus as WideT))
            .map(|v| v as ValueT)
            .collect();

        assert_eq!(shoup_res, normal_res);
    }
}
