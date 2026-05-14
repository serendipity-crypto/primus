//! Slice-level (bulk) Barrett operations.
//!
//! The scalar helpers in this module loop over a slice using the existing
//! scalar `Reduce*` impls on [`BarrettModulus`]. They are also used as
//! the tail of the SIMD kernels in [`super::simd`] when the
//! `nightly` + `simd` feature combo is enabled.

use primus_integer::UnsignedInteger;
use primus_reduce::{
    LazyReduceMul, LazyReduceMulAdd, LazyReduceMulAddSlice, LazyReduceMulSlice, ReduceAdd,
    ReduceAddAssign, ReduceAddSlice, ReduceMul, ReduceMulAdd, ReduceMulAddSlice, ReduceMulSlice,
    ReduceNeg, ReduceNegAssign, ReduceNegSlice, ReduceOnce, ReduceOnceAssign, ReduceOnceSlice,
    ReduceSub, ReduceSubAssign, ReduceSubSlice,
};

use super::BarrettModulus;

// ---------------------------------------------------------------------------
// Scalar slice helpers (also reused as the SIMD tail).
// ---------------------------------------------------------------------------

#[inline]
pub(super) fn scalar_reduce_once_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    values: &mut [T],
) {
    values
        .iter_mut()
        .for_each(|v| modulus.reduce_once_assign(v));
}

#[inline]
pub(super) fn scalar_reduce_once_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    input: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(input.len(), output.len());
    input
        .iter()
        .zip(output)
        .for_each(|(&v, o)| *o = modulus.reduce_once(v));
}

#[inline]
pub(super) fn scalar_reduce_neg_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    values: &mut [T],
) {
    values.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
}

#[inline]
pub(super) fn scalar_reduce_neg_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    input: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(input.len(), output.len());
    input
        .iter()
        .zip(output)
        .for_each(|(&v, o)| *o = modulus.reduce_neg(v));
}

#[inline]
pub(super) fn scalar_reduce_add_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) {
    debug_assert_eq!(a.len(), b.len());
    a.iter_mut()
        .zip(b)
        .for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
}

#[inline]
pub(super) fn scalar_reduce_add_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(output)
        .for_each(|((&a, &b), o)| *o = modulus.reduce_add(a, b));
}

#[inline]
pub(super) fn scalar_reduce_sub_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) {
    debug_assert_eq!(a.len(), b.len());
    a.iter_mut()
        .zip(b)
        .for_each(|(a, &b)| modulus.reduce_sub_assign(a, b));
}

#[inline]
pub(super) fn scalar_reduce_sub_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(output)
        .for_each(|((&a, &b), o)| *o = modulus.reduce_sub(a, b));
}

#[inline]
pub(super) fn scalar_reduce_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) {
    debug_assert_eq!(a.len(), b.len());
    a.iter_mut()
        .zip(b)
        .for_each(|(a, &b)| *a = modulus.reduce_mul(*a, b));
}

#[inline]
pub(super) fn scalar_reduce_mul_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(output)
        .for_each(|((&a, &b), o)| *o = modulus.reduce_mul(a, b));
}

#[inline]
pub(super) fn scalar_lazy_reduce_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) {
    debug_assert_eq!(a.len(), b.len());
    a.iter_mut()
        .zip(b)
        .for_each(|(a, &b)| *a = modulus.lazy_reduce_mul(*a, b));
}

#[inline]
pub(super) fn scalar_lazy_reduce_mul_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(output)
        .for_each(|((&a, &b), o)| *o = modulus.lazy_reduce_mul(a, b));
}

#[inline]
pub(super) fn scalar_reduce_add_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    acc.iter_mut()
        .zip(a)
        .zip(b)
        .for_each(|((acc, &a), &b)| *acc = modulus.reduce_mul_add(a, b, *acc));
}

#[inline]
pub(super) fn scalar_reduce_sub_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    acc.iter_mut().zip(a).zip(b).for_each(|((acc, &a), &b)| {
        let prod = modulus.reduce_mul(a, b);
        modulus.reduce_sub_assign(acc, prod);
    });
}

#[inline]
pub(super) fn scalar_reduce_mul_add_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    c: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(c)
        .zip(output)
        .for_each(|(((&a, &b), &c), o)| *o = modulus.reduce_mul_add(a, b, c));
}

#[inline]
pub(super) fn scalar_reduce_scalar_mul_add_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    scalar: T,
    b: &[T],
    c: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    b.iter()
        .zip(c)
        .zip(output)
        .for_each(|((&b, &c), o)| *o = modulus.reduce_mul_add(scalar, b, c));
}

#[inline]
pub(super) fn scalar_lazy_reduce_add_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    acc.iter_mut()
        .zip(a)
        .zip(b)
        .for_each(|((acc, &a), &b)| *acc = modulus.lazy_reduce_mul_add(a, b, *acc));
}

#[inline]
pub(super) fn scalar_lazy_reduce_sub_mul_slice_assign<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) {
    // Lazy here means: result is allowed in [0, 2*modulus). We compute
    // `acc - a*b` by rewriting it as `acc + (modulus - a*b mod modulus)`
    // so the intermediate stays unsigned, then keep the canonical-add
    // path; tightening to [0, 2m) is only profitable when no further
    // borrow is needed. We just route to canonical for now.
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    acc.iter_mut().zip(a).zip(b).for_each(|((acc, &a), &b)| {
        let prod = modulus.reduce_mul(a, b);
        modulus.reduce_sub_assign(acc, prod);
    });
}

#[inline]
pub(super) fn scalar_lazy_reduce_mul_add_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    c: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    a.iter()
        .zip(b)
        .zip(c)
        .zip(output)
        .for_each(|(((&a, &b), &c), o)| *o = modulus.lazy_reduce_mul_add(a, b, c));
}

#[inline]
pub(super) fn scalar_lazy_reduce_scalar_mul_add_slice_to<T: UnsignedInteger>(
    modulus: BarrettModulus<T>,
    scalar: T,
    b: &[T],
    c: &[T],
    output: &mut [T],
) {
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    b.iter()
        .zip(c)
        .zip(output)
        .for_each(|((&b, &c), o)| *o = modulus.lazy_reduce_mul_add(scalar, b, c));
}

// ---------------------------------------------------------------------------
// Trait impls — wires `BarrettModulus<T>` to either the scalar helpers or
// the SIMD kernels in `super::simd`, selected per type via the macros below.
// ---------------------------------------------------------------------------

macro_rules! impl_barrett_slice_scalar {
    ($($t:ty),* $(,)?) => {
        $(
            impl ReduceOnceSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_once_slice_assign(self, values: &mut [$t]) {
                    scalar_reduce_once_slice_assign(self, values)
                }
                #[inline]
                fn reduce_once_slice_to(self, input: &[$t], output: &mut [$t]) {
                    scalar_reduce_once_slice_to(self, input, output)
                }
            }

            impl ReduceNegSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_neg_slice_assign(self, values: &mut [$t]) {
                    scalar_reduce_neg_slice_assign(self, values)
                }
                #[inline]
                fn reduce_neg_slice_to(self, input: &[$t], output: &mut [$t]) {
                    scalar_reduce_neg_slice_to(self, input, output)
                }
            }

            impl ReduceAddSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_add_slice_assign(self, a: &mut [$t], b: &[$t]) {
                    scalar_reduce_add_slice_assign(self, a, b)
                }
                #[inline]
                fn reduce_add_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                    scalar_reduce_add_slice_to(self, a, b, output)
                }
            }

            impl ReduceSubSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_sub_slice_assign(self, a: &mut [$t], b: &[$t]) {
                    scalar_reduce_sub_slice_assign(self, a, b)
                }
                #[inline]
                fn reduce_sub_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                    scalar_reduce_sub_slice_to(self, a, b, output)
                }
            }

            impl ReduceMulSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_mul_slice_assign(self, a: &mut [$t], b: &[$t]) {
                    scalar_reduce_mul_slice_assign(self, a, b)
                }
                #[inline]
                fn reduce_mul_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                    scalar_reduce_mul_slice_to(self, a, b, output)
                }
            }

            impl LazyReduceMulSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn lazy_reduce_mul_slice_assign(self, a: &mut [$t], b: &[$t]) {
                    scalar_lazy_reduce_mul_slice_assign(self, a, b)
                }
                #[inline]
                fn lazy_reduce_mul_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                    scalar_lazy_reduce_mul_slice_to(self, a, b, output)
                }
            }

            impl ReduceMulAddSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn reduce_add_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                    scalar_reduce_add_mul_slice_assign(self, acc, a, b)
                }
                #[inline]
                fn reduce_sub_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                    scalar_reduce_sub_mul_slice_assign(self, acc, a, b)
                }
                #[inline]
                fn reduce_mul_add_slice_to(
                    self,
                    a: &[$t],
                    b: &[$t],
                    c: &[$t],
                    output: &mut [$t],
                ) {
                    scalar_reduce_mul_add_slice_to(self, a, b, c, output)
                }
                #[inline]
                fn reduce_scalar_mul_add_slice_to(
                    self,
                    scalar: $t,
                    b: &[$t],
                    c: &[$t],
                    output: &mut [$t],
                ) {
                    scalar_reduce_scalar_mul_add_slice_to(self, scalar, b, c, output)
                }
            }

            impl LazyReduceMulAddSlice<$t> for BarrettModulus<$t> {
                #[inline]
                fn lazy_reduce_add_mul_slice_assign(
                    self,
                    acc: &mut [$t],
                    a: &[$t],
                    b: &[$t],
                ) {
                    scalar_lazy_reduce_add_mul_slice_assign(self, acc, a, b)
                }
                #[inline]
                fn lazy_reduce_sub_mul_slice_assign(
                    self,
                    acc: &mut [$t],
                    a: &[$t],
                    b: &[$t],
                ) {
                    scalar_lazy_reduce_sub_mul_slice_assign(self, acc, a, b)
                }
                #[inline]
                fn lazy_reduce_mul_add_slice_to(
                    self,
                    a: &[$t],
                    b: &[$t],
                    c: &[$t],
                    output: &mut [$t],
                ) {
                    scalar_lazy_reduce_mul_add_slice_to(self, a, b, c, output)
                }
                #[inline]
                fn lazy_reduce_scalar_mul_add_slice_to(
                    self,
                    scalar: $t,
                    b: &[$t],
                    c: &[$t],
                    output: &mut [$t],
                ) {
                    scalar_lazy_reduce_scalar_mul_add_slice_to(self, scalar, b, c, output)
                }
            }
        )*
    };
}

// `u128` always falls back to scalar — `SimdBarrettModulus<u128, _>` is
// not viable on most targets.
impl_barrett_slice_scalar!(u128);

// When SIMD is off, every primitive width falls back to scalar.
#[cfg(not(all(feature = "nightly", feature = "simd")))]
impl_barrett_slice_scalar!(u8, u16, u32, u64, usize);

#[cfg(all(feature = "nightly", feature = "simd"))]
macro_rules! impl_barrett_slice_simd {
    ($t:ty, $lanes:expr) => {
        impl ReduceOnceSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_once_slice_assign(self, values: &mut [$t]) {
                super::simd::reduce_once_slice_assign::<$t, { $lanes }>(self, values)
            }
            #[inline]
            fn reduce_once_slice_to(self, input: &[$t], output: &mut [$t]) {
                super::simd::reduce_once_slice_to::<$t, { $lanes }>(self, input, output)
            }
        }

        impl ReduceNegSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_neg_slice_assign(self, values: &mut [$t]) {
                super::simd::reduce_neg_slice_assign::<$t, { $lanes }>(self, values)
            }
            #[inline]
            fn reduce_neg_slice_to(self, input: &[$t], output: &mut [$t]) {
                super::simd::reduce_neg_slice_to::<$t, { $lanes }>(self, input, output)
            }
        }

        impl ReduceAddSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_add_slice_assign(self, a: &mut [$t], b: &[$t]) {
                super::simd::reduce_add_slice_assign::<$t, { $lanes }>(self, a, b)
            }
            #[inline]
            fn reduce_add_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                super::simd::reduce_add_slice_to::<$t, { $lanes }>(self, a, b, output)
            }
        }

        impl ReduceSubSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_sub_slice_assign(self, a: &mut [$t], b: &[$t]) {
                super::simd::reduce_sub_slice_assign::<$t, { $lanes }>(self, a, b)
            }
            #[inline]
            fn reduce_sub_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                super::simd::reduce_sub_slice_to::<$t, { $lanes }>(self, a, b, output)
            }
        }

        impl ReduceMulSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_mul_slice_assign(self, a: &mut [$t], b: &[$t]) {
                super::simd::reduce_mul_slice_assign::<$t, { $lanes }>(self, a, b)
            }
            #[inline]
            fn reduce_mul_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                super::simd::reduce_mul_slice_to::<$t, { $lanes }>(self, a, b, output)
            }
        }

        impl LazyReduceMulSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn lazy_reduce_mul_slice_assign(self, a: &mut [$t], b: &[$t]) {
                super::simd::lazy_reduce_mul_slice_assign::<$t, { $lanes }>(self, a, b)
            }
            #[inline]
            fn lazy_reduce_mul_slice_to(self, a: &[$t], b: &[$t], output: &mut [$t]) {
                super::simd::lazy_reduce_mul_slice_to::<$t, { $lanes }>(self, a, b, output)
            }
        }

        impl ReduceMulAddSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn reduce_add_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                super::simd::reduce_add_mul_slice_assign::<$t, { $lanes }>(self, acc, a, b)
            }
            #[inline]
            fn reduce_sub_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                super::simd::reduce_sub_mul_slice_assign::<$t, { $lanes }>(self, acc, a, b)
            }
            #[inline]
            fn reduce_mul_add_slice_to(self, a: &[$t], b: &[$t], c: &[$t], output: &mut [$t]) {
                super::simd::reduce_mul_add_slice_to::<$t, { $lanes }>(self, a, b, c, output)
            }
            #[inline]
            fn reduce_scalar_mul_add_slice_to(
                self,
                scalar: $t,
                b: &[$t],
                c: &[$t],
                output: &mut [$t],
            ) {
                super::simd::reduce_scalar_mul_add_slice_to::<$t, { $lanes }>(
                    self, scalar, b, c, output,
                )
            }
        }

        impl LazyReduceMulAddSlice<$t> for BarrettModulus<$t> {
            #[inline]
            fn lazy_reduce_add_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                super::simd::lazy_reduce_add_mul_slice_assign::<$t, { $lanes }>(self, acc, a, b)
            }
            #[inline]
            fn lazy_reduce_sub_mul_slice_assign(self, acc: &mut [$t], a: &[$t], b: &[$t]) {
                super::simd::lazy_reduce_sub_mul_slice_assign::<$t, { $lanes }>(self, acc, a, b)
            }
            #[inline]
            fn lazy_reduce_mul_add_slice_to(self, a: &[$t], b: &[$t], c: &[$t], output: &mut [$t]) {
                super::simd::lazy_reduce_mul_add_slice_to::<$t, { $lanes }>(self, a, b, c, output)
            }
            #[inline]
            fn lazy_reduce_scalar_mul_add_slice_to(
                self,
                scalar: $t,
                b: &[$t],
                c: &[$t],
                output: &mut [$t],
            ) {
                super::simd::lazy_reduce_scalar_mul_add_slice_to::<$t, { $lanes }>(
                    self, scalar, b, c, output,
                )
            }
        }
    };
}

#[cfg(all(feature = "nightly", feature = "simd"))]
impl_barrett_slice_simd!(u8, primus_integer::default_lanes::VECTOR_BITS / 8);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_barrett_slice_simd!(u16, primus_integer::default_lanes::VECTOR_BITS / 16);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_barrett_slice_simd!(u32, primus_integer::default_lanes::VECTOR_BITS / 32);
#[cfg(all(feature = "nightly", feature = "simd"))]
impl_barrett_slice_simd!(u64, primus_integer::default_lanes::VECTOR_BITS / 64);

#[cfg(all(feature = "nightly", feature = "simd", target_pointer_width = "64"))]
impl_barrett_slice_simd!(usize, primus_integer::default_lanes::VECTOR_BITS / 64);
#[cfg(all(feature = "nightly", feature = "simd", target_pointer_width = "32"))]
impl_barrett_slice_simd!(usize, primus_integer::default_lanes::VECTOR_BITS / 32);

#[cfg(test)]
mod tests {
    use rand::distr::{Distribution, Uniform};

    use super::*;

    type V = u32;
    type W = u64;
    const MODULUS: V = 132_120_577;

    // Pick a length that is intentionally not a multiple of any reasonable
    // lane count, so the SIMD tail is exercised.
    const LEN: usize = 67;

    fn rand_slice(len: usize) -> Vec<V> {
        let mut rng = rand::rng();
        let distr = Uniform::new(0, MODULUS).unwrap();
        distr.sample_iter(&mut rng).take(len).collect()
    }

    fn mmod(x: V) -> V {
        x % MODULUS
    }

    fn mul_mod(a: V, b: V) -> V {
        ((a as W * b as W) % MODULUS as W) as V
    }

    #[test]
    fn once_slice_assign_and_to() {
        let m = BarrettModulus::<V>::new(MODULUS);
        // Inputs in [0, 2m).
        let data: Vec<V> = (0..LEN)
            .map(|i| ((i as V).wrapping_mul(7919)) % (MODULUS * 2 - 1))
            .collect();
        let expected: Vec<V> = data.iter().map(|&v| mmod(v)).collect();

        let mut assign = data.clone();
        m.reduce_once_slice_assign(&mut assign);
        assert_eq!(assign, expected);

        let mut to = vec![0; LEN];
        m.reduce_once_slice_to(&data, &mut to);
        assert_eq!(to, expected);
    }

    #[test]
    fn neg_slice_assign_and_to() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let data = rand_slice(LEN);
        let expected: Vec<V> = data
            .iter()
            .map(|&v| if v == 0 { 0 } else { MODULUS - v })
            .collect();

        let mut assign = data.clone();
        m.reduce_neg_slice_assign(&mut assign);
        assert_eq!(assign, expected);

        let mut to = vec![0; LEN];
        m.reduce_neg_slice_to(&data, &mut to);
        assert_eq!(to, expected);
    }

    #[test]
    fn add_slice_assign_and_to() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = a.iter().zip(&b).map(|(&a, &b)| mmod(a + b)).collect();

        let mut assign = a.clone();
        m.reduce_add_slice_assign(&mut assign, &b);
        assert_eq!(assign, expected);

        let mut to = vec![0; LEN];
        m.reduce_add_slice_to(&a, &b, &mut to);
        assert_eq!(to, expected);
    }

    #[test]
    fn sub_slice_assign_and_to() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = a
            .iter()
            .zip(&b)
            .map(|(&a, &b)| if a >= b { a - b } else { a + MODULUS - b })
            .collect();

        let mut assign = a.clone();
        m.reduce_sub_slice_assign(&mut assign, &b);
        assert_eq!(assign, expected);

        let mut to = vec![0; LEN];
        m.reduce_sub_slice_to(&a, &b, &mut to);
        assert_eq!(to, expected);
    }

    #[test]
    fn mul_slice_assign_and_to() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = a.iter().zip(&b).map(|(&a, &b)| mul_mod(a, b)).collect();

        let mut assign = a.clone();
        m.reduce_mul_slice_assign(&mut assign, &b);
        assert_eq!(assign, expected);

        let mut to = vec![0; LEN];
        m.reduce_mul_slice_to(&a, &b, &mut to);
        assert_eq!(to, expected);
    }

    #[test]
    fn lazy_mul_slice_matches_canonical_mod_m() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = a.iter().zip(&b).map(|(&a, &b)| mul_mod(a, b)).collect();

        let mut to = vec![0; LEN];
        m.lazy_reduce_mul_slice_to(&a, &b, &mut to);
        // Lazy result is in [0, 2m); fold and compare to canonical.
        for v in to.iter_mut() {
            if *v >= MODULUS {
                *v -= MODULUS;
            }
        }
        assert_eq!(to, expected);
    }

    #[test]
    fn add_mul_slice_assign_matches_scalar() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let acc = rand_slice(LEN);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = acc
            .iter()
            .zip(&a)
            .zip(&b)
            .map(|((&acc, &a), &b)| mmod(acc + mul_mod(a, b)))
            .collect();

        let mut acc_v = acc.clone();
        m.reduce_add_mul_slice_assign(&mut acc_v, &a, &b);
        assert_eq!(acc_v, expected);
    }

    #[test]
    fn sub_mul_slice_assign_matches_scalar() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let acc = rand_slice(LEN);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = acc
            .iter()
            .zip(&a)
            .zip(&b)
            .map(|((&acc, &a), &b)| {
                let prod = mul_mod(a, b);
                if acc >= prod {
                    acc - prod
                } else {
                    acc + MODULUS - prod
                }
            })
            .collect();

        let mut acc_v = acc.clone();
        m.reduce_sub_mul_slice_assign(&mut acc_v, &a, &b);
        assert_eq!(acc_v, expected);
    }

    #[test]
    fn mul_add_slice_to_matches_scalar() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let c = rand_slice(LEN);
        let expected: Vec<V> = a
            .iter()
            .zip(&b)
            .zip(&c)
            .map(|((&a, &b), &c)| mmod(mul_mod(a, b) + c))
            .collect();

        let mut out = vec![0; LEN];
        m.reduce_mul_add_slice_to(&a, &b, &c, &mut out);
        assert_eq!(out, expected);
    }

    #[test]
    fn scalar_mul_add_slice_to_matches_scalar() {
        let mut rng = rand::rng();
        let m = BarrettModulus::<V>::new(MODULUS);
        let scalar = Uniform::new(0, MODULUS).unwrap().sample(&mut rng);
        let b = rand_slice(LEN);
        let c = rand_slice(LEN);
        let expected: Vec<V> = b
            .iter()
            .zip(&c)
            .map(|(&b, &c)| mmod(mul_mod(scalar, b) + c))
            .collect();

        let mut out = vec![0; LEN];
        m.reduce_scalar_mul_add_slice_to(scalar, &b, &c, &mut out);
        assert_eq!(out, expected);
    }

    #[test]
    fn lazy_add_mul_slice_matches_canonical_mod_m() {
        let m = BarrettModulus::<V>::new(MODULUS);
        let acc = rand_slice(LEN);
        let a = rand_slice(LEN);
        let b = rand_slice(LEN);
        let expected: Vec<V> = acc
            .iter()
            .zip(&a)
            .zip(&b)
            .map(|((&acc, &a), &b)| mmod(acc + mul_mod(a, b)))
            .collect();

        let mut acc_v = acc.clone();
        m.lazy_reduce_add_mul_slice_assign(&mut acc_v, &a, &b);
        for v in acc_v.iter_mut() {
            if *v >= MODULUS {
                *v -= MODULUS;
            }
        }
        assert_eq!(acc_v, expected);
    }

    /// Exercise a length not divisible by the default lane count so the
    /// SIMD path's scalar tail is hit.
    #[test]
    fn tail_handling() {
        let m = BarrettModulus::<V>::new(MODULUS);
        for &len in &[1usize, 3, 7, 8, 9, 15, 16, 17, 31, 33, 64, 65] {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let expected: Vec<V> = a.iter().zip(&b).map(|(&a, &b)| mul_mod(a, b)).collect();

            let mut out = vec![0; len];
            m.reduce_mul_slice_to(&a, &b, &mut out);
            assert_eq!(out, expected, "mismatch at len={len}");
        }
    }
}
