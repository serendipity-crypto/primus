use primus_integer::{DivRemScalar, UnsignedInteger};

use crate::{FactorMul, FactorSliceOps, LazyFactorMul, LazyFactorSliceOps};

#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::SimdShoupFactor;

/// Lower-level SIMD slice kernels with an explicit lane count, for callers
/// that want to override the default vector width.
///
/// The default trait impls in [`ShoupFactorSliceOps`] pick a lane count at
/// compile time based on the target CPU's SIMD width
/// (see [`primus_integer::default_lanes`]). Reach for this module only when
/// you have measured a different `N` that performs better on your workload.
#[cfg(all(feature = "nightly", feature = "simd"))]
pub mod simd_kernel {
    pub use super::simd::{
        add_factor_mul_slice_assign, factor_mul_add_slice_to, factor_mul_slice_assign,
        factor_mul_slice_to, lazy_factor_mul_slice_assign, lazy_factor_mul_slice_to,
        sub_factor_mul_slice_assign,
    };
}

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

#[inline]
fn reduce_add<T: UnsignedInteger>(a: T, b: T, modulus: T) -> T {
    let sum = a + b;
    if sum >= modulus { sum - modulus } else { sum }
}

#[inline]
pub(super) fn scalar_lazy_factor_mul_slice_assign<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    values: &mut [T],
    modulus: T,
) {
    values
        .iter_mut()
        .for_each(|value| *value = factor.lazy_factor_mul_modulo(*value, modulus));
}

#[inline]
pub(super) fn scalar_lazy_factor_mul_slice_to<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    input: &[T],
    output: &mut [T],
    modulus: T,
) {
    assert_eq!(input.len(), output.len());
    input
        .iter()
        .zip(output)
        .for_each(|(&value, output)| *output = factor.lazy_factor_mul_modulo(value, modulus));
}

#[inline]
pub(super) fn scalar_factor_mul_slice_assign<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    values: &mut [T],
    modulus: T,
) {
    values
        .iter_mut()
        .for_each(|value| *value = factor.factor_mul_modulo(*value, modulus));
}

#[inline]
pub(super) fn scalar_factor_mul_slice_to<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    input: &[T],
    output: &mut [T],
    modulus: T,
) {
    assert_eq!(input.len(), output.len());
    input
        .iter()
        .zip(output)
        .for_each(|(&value, output)| *output = factor.factor_mul_modulo(value, modulus));
}

#[inline]
pub(super) fn scalar_add_factor_mul_slice_assign<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    acc: &mut [T],
    rhs: &[T],
    modulus: T,
) {
    assert_eq!(acc.len(), rhs.len());
    acc.iter_mut().zip(rhs).for_each(|(acc, &rhs)| {
        *acc = reduce_add(*acc, factor.factor_mul_modulo(rhs, modulus), modulus);
    });
}

#[inline]
pub(super) fn scalar_sub_factor_mul_slice_assign<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    acc: &mut [T],
    rhs: &[T],
    modulus: T,
) {
    assert_eq!(acc.len(), rhs.len());
    acc.iter_mut().zip(rhs).for_each(|(acc, &rhs)| {
        let prod = factor.factor_mul_modulo(rhs, modulus);
        // `*acc - prod (mod modulus)`. Both are in `[0, modulus)`.
        *acc = if *acc >= prod {
            *acc - prod
        } else {
            acc.wrapping_add(modulus).wrapping_sub(prod)
        };
    });
}

#[inline]
pub(super) fn scalar_factor_mul_add_slice_to<T: UnsignedInteger>(
    factor: ShoupFactor<T>,
    b: &[T],
    c: &[T],
    output: &mut [T],
    modulus: T,
) {
    assert_eq!(b.len(), c.len());
    assert_eq!(b.len(), output.len());
    b.iter().zip(c).zip(output).for_each(|((&b, &c), o)| {
        let prod = factor.factor_mul_modulo(b, modulus);
        *o = reduce_add(prod, c, modulus);
    });
}

macro_rules! impl_shoup_factor_slice_ops_scalar {
    ($($t:ty),* $(,)?) => {
        $(
            impl LazyFactorSliceOps<$t> for ShoupFactor<$t> {
                #[inline]
                fn lazy_factor_mul_slice_assign(self, values: &mut [$t], modulus: $t) {
                    scalar_lazy_factor_mul_slice_assign(self, values, modulus);
                }

                #[inline]
                fn lazy_factor_mul_slice_to(self, input: &[$t], output: &mut [$t], modulus: $t) {
                    scalar_lazy_factor_mul_slice_to(self, input, output, modulus);
                }
            }

            impl FactorSliceOps<$t> for ShoupFactor<$t> {
                #[inline]
                fn factor_mul_slice_assign(self, values: &mut [$t], modulus: $t) {
                    scalar_factor_mul_slice_assign(self, values, modulus);
                }

                #[inline]
                fn factor_mul_slice_to(self, input: &[$t], output: &mut [$t], modulus: $t) {
                    scalar_factor_mul_slice_to(self, input, output, modulus);
                }

                #[inline]
                fn add_factor_mul_slice_assign(self, acc: &mut [$t], rhs: &[$t], modulus: $t) {
                    scalar_add_factor_mul_slice_assign(self, acc, rhs, modulus);
                }

                #[inline]
                fn sub_factor_mul_slice_assign(self, acc: &mut [$t], rhs: &[$t], modulus: $t) {
                    scalar_sub_factor_mul_slice_assign(self, acc, rhs, modulus);
                }

                #[inline]
                fn factor_mul_add_slice_to(
                    self,
                    b: &[$t],
                    c: &[$t],
                    output: &mut [$t],
                    modulus: $t,
                ) {
                    scalar_factor_mul_add_slice_to(self, b, c, output, modulus);
                }
            }
        )*
    };
}

#[cfg(all(feature = "nightly", feature = "simd"))]
macro_rules! impl_shoup_factor_slice_ops_simd {
    ($t:ty, $lanes:expr) => {
        impl LazyFactorSliceOps<$t> for ShoupFactor<$t> {
            #[inline]
            fn lazy_factor_mul_slice_assign(self, values: &mut [$t], modulus: $t) {
                simd::lazy_factor_mul_slice_assign::<$t, { $lanes }>(self, values, modulus);
            }

            #[inline]
            fn lazy_factor_mul_slice_to(self, input: &[$t], output: &mut [$t], modulus: $t) {
                simd::lazy_factor_mul_slice_to::<$t, { $lanes }>(self, input, output, modulus);
            }
        }

        impl FactorSliceOps<$t> for ShoupFactor<$t> {
            #[inline]
            fn factor_mul_slice_assign(self, values: &mut [$t], modulus: $t) {
                simd::factor_mul_slice_assign::<$t, { $lanes }>(self, values, modulus);
            }

            #[inline]
            fn factor_mul_slice_to(self, input: &[$t], output: &mut [$t], modulus: $t) {
                simd::factor_mul_slice_to::<$t, { $lanes }>(self, input, output, modulus);
            }

            #[inline]
            fn add_factor_mul_slice_assign(self, acc: &mut [$t], rhs: &[$t], modulus: $t) {
                simd::add_factor_mul_slice_assign::<$t, { $lanes }>(self, acc, rhs, modulus);
            }

            #[inline]
            fn sub_factor_mul_slice_assign(self, acc: &mut [$t], rhs: &[$t], modulus: $t) {
                simd::sub_factor_mul_slice_assign::<$t, { $lanes }>(self, acc, rhs, modulus);
            }

            #[inline]
            fn factor_mul_add_slice_to(
                self,
                b: &[$t],
                c: &[$t],
                output: &mut [$t],
                modulus: $t,
            ) {
                simd::factor_mul_add_slice_to::<$t, { $lanes }>(self, b, c, output, modulus);
            }
        }
    };
}

#[cfg(all(feature = "nightly", feature = "simd"))]
impl_shoup_factor_slice_ops_simd!(u8, primus_integer::default_lanes::VECTOR_BITS / 8);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_shoup_factor_slice_ops_simd!(u16, primus_integer::default_lanes::VECTOR_BITS / 16);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_shoup_factor_slice_ops_simd!(u32, primus_integer::default_lanes::VECTOR_BITS / 32);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_shoup_factor_slice_ops_simd!(u64, primus_integer::default_lanes::VECTOR_BITS / 64);

#[cfg(all(feature = "nightly", feature = "simd", target_pointer_width = "64"))]
impl_shoup_factor_slice_ops_simd!(usize, primus_integer::default_lanes::VECTOR_BITS / 64);
#[cfg(all(feature = "nightly", feature = "simd", target_pointer_width = "32"))]
impl_shoup_factor_slice_ops_simd!(usize, primus_integer::default_lanes::VECTOR_BITS / 32);

#[cfg(not(all(feature = "nightly", feature = "simd")))]
impl_shoup_factor_slice_ops_scalar!(u8, u16, u32, u64, usize);

impl_shoup_factor_slice_ops_scalar!(u128);

#[cfg(test)]
mod tests {
    use rand::{
        RngExt,
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

    #[test]
    fn test_shoup_slice_ops() {
        let mut rng = rand::rng();

        let modulus: ValueT = 132120577;
        let distr = Uniform::new(0, modulus).unwrap();

        let shoup = ShoupFactor::new(rng.sample(distr), modulus);
        let data: Vec<ValueT> = distr.sample_iter(&mut rng).take(N + 5).collect();

        let expected_lazy: Vec<ValueT> = data
            .iter()
            .map(|&v| shoup.lazy_factor_mul_modulo(v, modulus))
            .collect();
        let mut lazy_assign = data.clone();
        shoup.lazy_factor_mul_slice_assign(&mut lazy_assign, modulus);
        assert_eq!(lazy_assign, expected_lazy);

        let expected: Vec<ValueT> = data
            .iter()
            .map(|&v| shoup.factor_mul_modulo(v, modulus))
            .collect();

        let mut assign = data.clone();
        shoup.factor_mul_slice_assign(&mut assign, modulus);
        assert_eq!(assign, expected);

        let mut output = vec![ValueT::default(); data.len()];
        shoup.factor_mul_slice_to(&data, &mut output, modulus);
        assert_eq!(output, expected);

        let mut acc: Vec<ValueT> = distr.sample_iter(&mut rng).take(data.len()).collect();
        let expected_add: Vec<ValueT> = acc
            .iter()
            .zip(&data)
            .map(|(&acc, &value)| {
                ((acc as WideT + ((shoup.value as WideT * value as WideT) % modulus as WideT))
                    % modulus as WideT) as ValueT
            })
            .collect();

        shoup.add_factor_mul_slice_assign(&mut acc, &data, modulus);
        assert_eq!(acc, expected_add);

        // sub_factor_mul: acc -= factor * data (mod modulus)
        let acc_init: Vec<ValueT> = distr.sample_iter(&mut rng).take(data.len()).collect();
        let expected_sub: Vec<ValueT> = acc_init
            .iter()
            .zip(&data)
            .map(|(&acc, &value)| {
                let prod = ((shoup.value as WideT * value as WideT) % modulus as WideT) as ValueT;
                if acc >= prod {
                    acc - prod
                } else {
                    acc + modulus - prod
                }
            })
            .collect();
        let mut acc_sub = acc_init.clone();
        shoup.sub_factor_mul_slice_assign(&mut acc_sub, &data, modulus);
        assert_eq!(acc_sub, expected_sub);

        // factor_mul_add: out = factor * data + c (mod modulus)
        let c: Vec<ValueT> = distr.sample_iter(&mut rng).take(data.len()).collect();
        let expected_mul_add: Vec<ValueT> = data
            .iter()
            .zip(&c)
            .map(|(&v, &c)| {
                let prod = (shoup.value as WideT * v as WideT) % modulus as WideT;
                ((prod + c as WideT) % modulus as WideT) as ValueT
            })
            .collect();
        let mut out_mul_add = vec![ValueT::default(); data.len()];
        shoup.factor_mul_add_slice_to(&data, &c, &mut out_mul_add, modulus);
        assert_eq!(out_mul_add, expected_mul_add);
    }

    #[test]
    fn test_shoup_slice_ops_u128_scalar_backend() {
        let modulus = 1_125_899_906_826_241u128;
        let shoup = ShoupFactor::new(123_456_789u128, modulus);
        let data = [
            0u128,
            1,
            42,
            132_120_577,
            536_813_569,
            modulus - 1,
            modulus / 2,
        ];

        let expected: Vec<_> = data
            .iter()
            .map(|&value| shoup.factor_mul_modulo(value, modulus))
            .collect();
        let mut output = vec![0u128; data.len()];

        shoup.factor_mul_slice_to(&data, &mut output, modulus);
        assert_eq!(output, expected);
    }

    /// Verifies that the SIMD kernel works at a non-default lane count.
    #[cfg(all(feature = "nightly", feature = "simd"))]
    #[test]
    fn test_shoup_simd_kernel_custom_lanes() {
        let mut rng = rand::rng();

        let modulus: ValueT = 132120577;
        let distr = Uniform::new(0, modulus).unwrap();

        let shoup = ShoupFactor::new(rng.sample(distr), modulus);
        let data: Vec<ValueT> = distr.sample_iter(&mut rng).take(67).collect();

        let expected: Vec<ValueT> = data
            .iter()
            .map(|&v| shoup.factor_mul_modulo(v, modulus))
            .collect();

        // Force a different lane count than the compile-time default.
        let mut output = vec![ValueT::default(); data.len()];
        simd_kernel::factor_mul_slice_to::<ValueT, 16>(shoup, &data, &mut output, modulus);
        assert_eq!(output, expected);

        let mut output = vec![ValueT::default(); data.len()];
        simd_kernel::factor_mul_slice_to::<ValueT, 4>(shoup, &data, &mut output, modulus);
        assert_eq!(output, expected);
    }
}
