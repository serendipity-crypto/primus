use core::simd::{
    Simd,
    cmp::{SimdOrd, SimdPartialOrd},
};

use primus_integer::{DivRemScalar, SimdArray, SimdMaskArray, SimdUnsignedInteger, WideningMul};

use crate::{FactorMul, LazyFactorMul};

use super::ShoupFactor;

/// A number used for fast modular multiplication.
///
/// This is efficient if many operations are multiplied by
/// the same number and then reduced with the same modulus.
#[derive(Debug, Clone, Copy, Default)]
pub struct SimdShoupFactor<T: SimdUnsignedInteger, const N: usize> {
    /// value
    value: Simd<T, N>,
    /// quotient
    quotient: Simd<T, N>,
}

impl<T: SimdUnsignedInteger, const N: usize> From<ShoupFactor<T>> for SimdShoupFactor<T, N> {
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

impl<T: SimdUnsignedInteger, const N: usize> LazyFactorMul<Simd<T, N>> for SimdShoupFactor<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_factor_mul_modulo(self, b: Simd<T, N>, modulus: Simd<T, N>) -> Simd<T, N> {
        let hw = self.quotient.widening_mul_hw(b);
        self.value * b - (modulus * hw)
    }
}

impl<T: SimdUnsignedInteger, const N: usize> FactorMul<Simd<T, N>> for SimdShoupFactor<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn factor_mul_modulo(self, b: Simd<T, N>, modulus: Simd<T, N>) -> Simd<T, N> {
        let t = self.lazy_factor_mul_modulo(b, modulus);
        t.simd_ge(modulus).select(t - modulus, t)
    }
}

#[inline]
pub fn lazy_factor_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    values: &mut [T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (chunks, remainder) = values.as_chunks_mut::<N>();
    for chunk in chunks {
        let value = Simd::from_array(*chunk);
        *chunk = simd_factor
            .lazy_factor_mul_modulo(value, simd_modulus)
            .to_array();
    }

    super::scalar_lazy_factor_mul_slice_assign(factor, remainder, modulus);
}

#[inline]
pub fn lazy_factor_mul_slice_to<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    input: &[T],
    output: &mut [T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(input.len(), output.len());

    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (input_chunks, input_rem) = input.as_chunks::<N>();
    let (output_chunks, output_rem) = output.as_chunks_mut::<N>();
    for (input, output) in input_chunks.iter().zip(output_chunks) {
        let value = Simd::from_array(*input);
        *output = simd_factor
            .lazy_factor_mul_modulo(value, simd_modulus)
            .to_array();
    }

    super::scalar_lazy_factor_mul_slice_to(factor, input_rem, output_rem, modulus);
}

#[inline]
pub fn factor_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    values: &mut [T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (chunks, remainder) = values.as_chunks_mut::<N>();
    for chunk in chunks {
        let value = Simd::from_array(*chunk);
        *chunk = simd_factor
            .factor_mul_modulo(value, simd_modulus)
            .to_array();
    }

    super::scalar_factor_mul_slice_assign(factor, remainder, modulus);
}

#[inline]
pub fn factor_mul_slice_to<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    input: &[T],
    output: &mut [T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(input.len(), output.len());

    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (input_chunks, input_rem) = input.as_chunks::<N>();
    let (output_chunks, output_rem) = output.as_chunks_mut::<N>();
    for (input, output) in input_chunks.iter().zip(output_chunks) {
        let value = Simd::from_array(*input);
        *output = simd_factor
            .factor_mul_modulo(value, simd_modulus)
            .to_array();
    }

    super::scalar_factor_mul_slice_to(factor, input_rem, output_rem, modulus);
}

#[inline]
pub fn add_factor_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    acc: &mut [T],
    rhs: &[T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(acc.len(), rhs.len());

    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (rhs_chunks, rhs_rem) = rhs.as_chunks::<N>();
    for (acc, rhs) in acc_chunks.iter_mut().zip(rhs_chunks) {
        let acc_value = Simd::from_array(*acc);
        let product = simd_factor.factor_mul_modulo(Simd::from_array(*rhs), simd_modulus);
        let sum = acc_value + product;
        // `acc, product ∈ [0, modulus)` ⇒ `sum ∈ [0, 2*modulus)`. The
        // `min(sum, sum - m)` trick lowers to `vpminuq` on AVX-512.
        *acc = sum.simd_min(sum - simd_modulus).to_array();
    }

    super::scalar_add_factor_mul_slice_assign(factor, acc_rem, rhs_rem, modulus);
}

#[inline]
pub fn sub_factor_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    acc: &mut [T],
    rhs: &[T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(acc.len(), rhs.len());

    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (rhs_chunks, rhs_rem) = rhs.as_chunks::<N>();
    for (acc, rhs) in acc_chunks.iter_mut().zip(rhs_chunks) {
        let acc_value = Simd::from_array(*acc);
        let product = simd_factor.factor_mul_modulo(Simd::from_array(*rhs), simd_modulus);
        // Branchless `(acc - product) mod modulus`.
        let diff = acc_value - product;
        *acc = acc_value
            .simd_lt(product)
            .select(diff + simd_modulus, diff)
            .to_array();
    }

    super::scalar_sub_factor_mul_slice_assign(factor, acc_rem, rhs_rem, modulus);
}

#[inline]
pub fn factor_mul_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    factor: ShoupFactor<T>,
    b: &[T],
    c: &[T],
    output: &mut [T],
    modulus: T,
) where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(b.len(), c.len());
    assert_eq!(b.len(), output.len());

    let simd_factor = SimdShoupFactor::<T, N>::from(factor);
    let simd_modulus = Simd::splat(modulus);

    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (c_chunks, c_rem) = c.as_chunks::<N>();
    let (out_chunks, out_rem) = output.as_chunks_mut::<N>();
    for ((bc, cc), oc) in b_chunks.iter().zip(c_chunks).zip(out_chunks) {
        let bv = Simd::from_array(*bc);
        let cv = Simd::from_array(*cc);
        let product = simd_factor.factor_mul_modulo(bv, simd_modulus);
        // `product, c ∈ [0, modulus)` ⇒ `sum ∈ [0, 2*modulus)`.
        let sum = product + cv;
        *oc = sum.simd_min(sum - simd_modulus).to_array();
    }

    super::scalar_factor_mul_add_slice_to(factor, b_rem, c_rem, out_rem, modulus);
}
