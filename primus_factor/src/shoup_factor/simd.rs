use core::simd::{LaneCount, Simd, SupportedLaneCount, cmp::SimdPartialOrd};

use integer::{DivRemScalar, SimdArray, SimdMaskArray, SimdUnsignedInteger, WideningMul};

use crate::{FactorMul, LazyFactorMul};

use super::ShoupFactor;

/// A number used for fast modular multiplication.
///
/// This is efficient if many operations are multiplied by
/// the same number and then reduced with the same modulus.
#[derive(Debug, Clone, Copy, Default)]
pub struct SimdShoupFactor<T: SimdUnsignedInteger, const N: usize>
where
    LaneCount<N>: SupportedLaneCount,
{
    /// value
    value: Simd<T, N>,
    /// quotient
    quotient: Simd<T, N>,
}

impl<T: SimdUnsignedInteger, const N: usize> From<ShoupFactor<T>> for SimdShoupFactor<T, N>
where
    LaneCount<N>: SupportedLaneCount,
{
    #[inline]
    fn from(factor: ShoupFactor<T>) -> Self {
        Self {
            value: Simd::splat(factor.value()),
            quotient: Simd::splat(factor.quotient()),
        }
    }
}

impl<T: SimdUnsignedInteger, const N: usize> SimdShoupFactor<T, N>
where
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    /// Constructs a [`SimdShoupFactor<T, N>`].
    ///
    /// * `value` must be less than `modulus`.
    #[inline]
    pub fn new(value: Simd<T, N>, modulus: T) -> Self {
        debug_assert!(value.simd_lt(Simd::splat(modulus)).all());

        // Calculate the quotient of `value * 2^64 / modulus`.
        let mut quotient = [T::ZERO; N];
        let values = value.as_array();
        for i in 0..N {
            let mut temp = [T::ZERO; 2];
            DivRemScalar::div_rem_scalar(&[T::ZERO, values[i]], modulus, &mut temp);
            quotient[i] = temp[0];
        }

        Self {
            value,
            quotient: quotient.into(),
        }
    }

    /// Constructs a [`SimdShoupFactor<T, N>`].
    #[inline]
    pub fn with_quotient(value: Simd<T, N>, quotient: Simd<T, N>) -> Self {
        Self { value, quotient }
    }

    /// Resets the `modulus` of [`SimdShoupFactor<T, N>`].
    #[inline]
    pub fn set_modulus(&mut self, modulus: T) {
        debug_assert!(self.value.simd_lt(Simd::splat(modulus)).all());

        // Calculate the quotient of `value * 2^64 / modulus`.
        let mut quotient = [T::ZERO; N];
        let values = self.value.as_array();
        for i in 0..N {
            let mut temp = [T::ZERO; 2];
            DivRemScalar::div_rem_scalar(&[T::ZERO, values[i]], modulus, &mut temp);
            quotient[i] = temp[0];
        }

        self.quotient = quotient.into();
    }

    /// Resets the content of [`SimdShoupFactor<T, N>`].
    ///
    /// * `value` must be less than `modulus`.
    #[inline]
    pub fn set(&mut self, value: Simd<T, N>, modulus: T) {
        self.value = value;
        self.set_modulus(modulus);
    }

    /// Returns the value of this [`SimdShoupFactor<T, N>`].
    #[inline]
    pub fn value(self) -> Simd<T, N> {
        self.value
    }

    /// Returns the quotient of this [`SimdShoupFactor<T, N>`].
    #[inline]
    pub fn quotient(self) -> Simd<T, N> {
        self.quotient
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyFactorMul<Simd<T, N>, Simd<T, N>>
    for SimdShoupFactor<T, N>
where
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_factor_mul_modulo(self, b: Simd<T, N>, modulus: Simd<T, N>) -> Simd<T, N> {
        let hw = self.quotient.widening_mul_hw(b);
        self.value * b - (modulus * hw)
    }
}

impl<T: SimdUnsignedInteger, const N: usize> FactorMul<Simd<T, N>, Simd<T, N>>
    for SimdShoupFactor<T, N>
where
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn factor_mul_modulo(self, b: Simd<T, N>, modulus: Simd<T, N>) -> Simd<T, N> {
        let t = self.lazy_factor_mul_modulo(b, modulus);
        t.simd_ge(modulus).select(t - modulus, t)
    }
}
