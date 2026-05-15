use std::simd::{
    Simd,
    cmp::{SimdOrd, SimdPartialEq, SimdPartialOrd},
};

use primus_integer::{
    CarryingAdd, CarryingMul, SimdArray, SimdMaskArray, SimdUnsignedInteger, WideningMul,
};
use primus_reduce::{ReduceAdd, lazy_ops::*};

use super::BarrettModulus;

/// A modulus, using barrett reduction algorithm.
///
/// The struct stores the modulus number and some precomputed
/// data. Here, `b` = 2^T::BITS
///
/// It's efficient if many reductions are performed with a single modulus.
#[derive(Debug, Clone, Copy)]
pub struct SimdBarrettModulus<T: SimdUnsignedInteger, const N: usize>
where
    Simd<T, N>: SimdArray<T, N>,
{
    value: Simd<T, N>,
    ratio: [Simd<T, N>; 2],
}

impl<T: SimdUnsignedInteger, const N: usize> From<BarrettModulus<T>> for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn from(modulus: BarrettModulus<T>) -> Self {
        let ratio = modulus.ratio();
        Self {
            value: Simd::splat(modulus.value()),
            ratio: [Simd::splat(ratio[0]), Simd::splat(ratio[1])],
        }
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduce<Simd<T, N>> for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    type Output = Simd<T, N>;

    #[inline]
    fn lazy_reduce(self, value: Simd<T, N>) -> Self::Output {
        let tmp = value.widening_mul_hw(self.ratio[0]); // tmp1
        let q = value.carrying_mul_hw(self.ratio[1], tmp); // q₃

        // Step 2.
        value - (q * self.value) // r = r₁ - r₂
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduce<[Simd<T, N>; 2]>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    type Output = Simd<T, N>;

    #[inline]
    fn lazy_reduce(self, value: [Simd<T, N>; 2]) -> Self::Output {
        let ah = value[0].widening_mul_hw(self.ratio[0]);

        let b = value[0].carrying_mul(self.ratio[1], ah);
        let c = value[1].widening_mul(self.ratio[0]);

        let d = value[1] * self.ratio[1];

        let bch = b.1.carrying_add(c.1, b.0.overflowing_add(c.0).1).0;

        let q = d + bch;

        // Step 2.
        value[0] - (q * self.value)
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduce<(Simd<T, N>, Simd<T, N>)>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    type Output = Simd<T, N>;

    #[inline]
    fn lazy_reduce(self, value: (Simd<T, N>, Simd<T, N>)) -> Self::Output {
        let ah = value.0.widening_mul_hw(self.ratio[0]);

        let b = value.0.carrying_mul(self.ratio[1], ah);
        let c = value.1.widening_mul(self.ratio[0]);

        let d = value.1 * self.ratio[1];

        let bch = b.1.carrying_add(c.1, b.0.overflowing_add(c.0).1).0;

        let q = d + bch;

        // Step 2.
        value.0 - (q * self.value)
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceAssign<Simd<T, N>>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_reduce_assign(self, value: &mut Simd<T, N>) {
        *value = self.lazy_reduce(*value);
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceMul<Simd<T, N>> for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    type Output = Simd<T, N>;

    #[inline]
    fn lazy_reduce_mul(self, a: Simd<T, N>, b: Simd<T, N>) -> Self::Output {
        self.lazy_reduce(a.widening_mul(b))
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceMulAssign<Simd<T, N>>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_reduce_mul_assign(self, a: &mut Simd<T, N>, b: Simd<T, N>) {
        *a = self.lazy_reduce(a.widening_mul(b));
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceMulAdd<Simd<T, N>>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    type Output = Simd<T, N>;

    #[inline]
    fn lazy_reduce_mul_add(self, a: Simd<T, N>, b: Simd<T, N>, c: Simd<T, N>) -> Self::Output {
        self.lazy_reduce(a.carrying_mul(b, c))
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceMulAddAssign<Simd<T, N>>
    for SimdBarrettModulus<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_reduce_mul_add_assign(self, a: &mut Simd<T, N>, b: Simd<T, N>, c: Simd<T, N>) {
        *a = self.lazy_reduce(a.carrying_mul(b, c));
    }
}

// ===========================================================================
// SIMD slice kernels.
//
// Each kernel processes `N`-wide chunks with `SimdBarrettModulus` and
// delegates the remainder to the scalar slice helpers in `super::slice`.
// They are exposed as `pub fn` so callers can override the lane count
// (see `simd_kernel` in `mod.rs`).
// ===========================================================================

#[inline]
fn simd_reduce_once<T: SimdUnsignedInteger, const N: usize>(
    v: Simd<T, N>,
    m: Simd<T, N>,
) -> Simd<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    // `min(v, v - m)` trick: when `v < m`, `v - m` wraps to a huge value so
    // unsigned min picks `v`; when `v >= m`, `v - m` is the canonical form
    // and is smaller than `v`. Lowers to a single `vpminuq` on AVX-512
    // (vs. compare + blend + sub).
    v.simd_min(v - m)
}

#[inline]
fn simd_reduce_add<T: SimdUnsignedInteger, const N: usize>(
    a: Simd<T, N>,
    b: Simd<T, N>,
    m: Simd<T, N>,
) -> Simd<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    let sum = a + b;
    sum.simd_min(sum - m)
}

#[inline]
fn simd_reduce_sub<T: SimdUnsignedInteger, const N: usize>(
    a: Simd<T, N>,
    b: Simd<T, N>,
    m: Simd<T, N>,
) -> Simd<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    // Branchless: compute the wrapping difference first, then add `m` back
    // only on lanes where the original `a < b` borrowed. On AVX-512 this
    // tends to lower to `vpsubq` + `vpcmpltuq` + `vpaddq{k}` (masked add).
    let diff = a - b;
    a.simd_lt(b).select(diff + m, diff)
}

#[inline]
fn simd_reduce_neg<T: SimdUnsignedInteger, const N: usize>(
    v: Simd<T, N>,
    m: Simd<T, N>,
) -> Simd<T, N>
where
    Simd<T, N>: SimdArray<T, N>,
{
    let zero = Simd::splat(T::ZERO);
    v.simd_eq(zero).select(zero, m - v)
}

#[inline]
pub fn reduce_once_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    values: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    let m = Simd::splat(modulus.value());
    let (chunks, rem) = values.as_chunks_mut::<N>();
    for chunk in chunks {
        let v = Simd::from_array(*chunk);
        *chunk = simd_reduce_once(v, m).to_array();
    }
    super::slice::scalar_reduce_once_slice_assign(modulus, rem);
}

#[inline]
pub fn reduce_once_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    input: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(input.len(), output.len());
    let m = Simd::splat(modulus.value());
    let (in_chunks, in_rem) = input.as_chunks::<N>();
    let (out_chunks, out_rem) = output.as_chunks_mut::<N>();
    for (i, o) in in_chunks.iter().zip(out_chunks) {
        let v = Simd::from_array(*i);
        *o = simd_reduce_once(v, m).to_array();
    }
    super::slice::scalar_reduce_once_slice_to(modulus, in_rem, out_rem);
}

#[inline]
pub fn reduce_neg_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    values: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    let m = Simd::splat(modulus.value());
    let (chunks, rem) = values.as_chunks_mut::<N>();
    for chunk in chunks {
        let v = Simd::from_array(*chunk);
        *chunk = simd_reduce_neg(v, m).to_array();
    }
    super::slice::scalar_reduce_neg_slice_assign(modulus, rem);
}

#[inline]
pub fn reduce_neg_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    input: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(input.len(), output.len());
    let m = Simd::splat(modulus.value());
    let (in_chunks, in_rem) = input.as_chunks::<N>();
    let (out_chunks, out_rem) = output.as_chunks_mut::<N>();
    for (i, o) in in_chunks.iter().zip(out_chunks) {
        let v = Simd::from_array(*i);
        *o = simd_reduce_neg(v, m).to_array();
    }
    super::slice::scalar_reduce_neg_slice_to(modulus, in_rem, out_rem);
}

#[inline]
pub fn reduce_add_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks_mut::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for (ac, bc) in a_chunks.iter_mut().zip(b_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *ac = simd_reduce_add(av, bv, m).to_array();
    }
    super::slice::scalar_reduce_add_slice_assign(modulus, a_rem, b_rem);
}

#[inline]
pub fn reduce_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((ac, bc), oc) in a_chunks.iter().zip(b_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *oc = simd_reduce_add(av, bv, m).to_array();
    }
    super::slice::scalar_reduce_add_slice_to(modulus, a_rem, b_rem, o_rem);
}

#[inline]
pub fn reduce_sub_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks_mut::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for (ac, bc) in a_chunks.iter_mut().zip(b_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *ac = simd_reduce_sub(av, bv, m).to_array();
    }
    super::slice::scalar_reduce_sub_slice_assign(modulus, a_rem, b_rem);
}

#[inline]
pub fn reduce_sub_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((ac, bc), oc) in a_chunks.iter().zip(b_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *oc = simd_reduce_sub(av, bv, m).to_array();
    }
    super::slice::scalar_reduce_sub_slice_to(modulus, a_rem, b_rem, o_rem);
}

#[inline]
pub fn lazy_reduce_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let (a_chunks, a_rem) = a.as_chunks_mut::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for (ac, bc) in a_chunks.iter_mut().zip(b_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *ac = sm.lazy_reduce_mul(av, bv).to_array();
    }
    super::slice::scalar_lazy_reduce_mul_slice_assign(modulus, a_rem, b_rem);
}

#[inline]
pub fn lazy_reduce_mul_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((ac, bc), oc) in a_chunks.iter().zip(b_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *oc = sm.lazy_reduce_mul(av, bv).to_array();
    }
    super::slice::scalar_lazy_reduce_mul_slice_to(modulus, a_rem, b_rem, o_rem);
}

#[inline]
pub fn reduce_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &mut [T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks_mut::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for (ac, bc) in a_chunks.iter_mut().zip(b_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let lazy = sm.lazy_reduce_mul(av, bv);
        *ac = simd_reduce_once(lazy, m).to_array();
    }
    super::slice::scalar_reduce_mul_slice_assign(modulus, a_rem, b_rem);
}

#[inline]
pub fn reduce_mul_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((ac, bc), oc) in a_chunks.iter().zip(b_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let lazy = sm.lazy_reduce_mul(av, bv);
        *oc = simd_reduce_once(lazy, m).to_array();
    }
    super::slice::scalar_reduce_mul_slice_to(modulus, a_rem, b_rem, o_rem);
}

#[inline]
pub fn reduce_add_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for ((accc, ac), bc) in acc_chunks.iter_mut().zip(a_chunks).zip(b_chunks) {
        let accv = Simd::from_array(*accc);
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        // (a*b + acc) mod 2m is in [0, 2m); fold once to canonical.
        let lazy = sm.lazy_reduce_mul_add(av, bv, accv);
        *accc = simd_reduce_once(lazy, m).to_array();
    }
    super::slice::scalar_reduce_add_mul_slice_assign(modulus, acc_rem, a_rem, b_rem);
}

#[inline]
pub fn reduce_sub_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for ((accc, ac), bc) in acc_chunks.iter_mut().zip(a_chunks).zip(b_chunks) {
        let accv = Simd::from_array(*accc);
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let prod_lazy = sm.lazy_reduce_mul(av, bv);
        let prod = simd_reduce_once(prod_lazy, m);
        *accc = simd_reduce_sub(accv, prod, m).to_array();
    }
    super::slice::scalar_reduce_sub_mul_slice_assign(modulus, acc_rem, a_rem, b_rem);
}

#[inline]
pub fn reduce_mul_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    c: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (c_chunks, c_rem) = c.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for (((ac, bc), cc), oc) in a_chunks.iter().zip(b_chunks).zip(c_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let cv = Simd::from_array(*cc);
        let lazy = sm.lazy_reduce_mul_add(av, bv, cv);
        *oc = simd_reduce_once(lazy, m).to_array();
    }
    super::slice::scalar_reduce_mul_add_slice_to(modulus, a_rem, b_rem, c_rem, o_rem);
}

#[inline]
pub fn reduce_scalar_mul_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    scalar: T,
    b: &[T],
    c: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let sv = Simd::splat(scalar);
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (c_chunks, c_rem) = c.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((bc, cc), oc) in b_chunks.iter().zip(c_chunks).zip(o_chunks) {
        let bv = Simd::from_array(*bc);
        let cv = Simd::from_array(*cc);
        let lazy = sm.lazy_reduce_mul_add(sv, bv, cv);
        *oc = simd_reduce_once(lazy, m).to_array();
    }
    super::slice::scalar_reduce_scalar_mul_add_slice_to(modulus, scalar, b_rem, c_rem, o_rem);
}

#[inline]
pub fn lazy_reduce_add_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for ((accc, ac), bc) in acc_chunks.iter_mut().zip(a_chunks).zip(b_chunks) {
        let accv = Simd::from_array(*accc);
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        *accc = sm.lazy_reduce_mul_add(av, bv, accv).to_array();
    }
    super::slice::scalar_lazy_reduce_add_mul_slice_assign(modulus, acc_rem, a_rem, b_rem);
}

#[inline]
pub fn lazy_reduce_sub_mul_slice_assign<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    acc: &mut [T],
    a: &[T],
    b: &[T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    // True lazy form: skip the `reduce_once` on the product.
    //
    // `prod_lazy ∈ [0, 2m)`, `acc ∈ [0, m)`. We want
    // `(acc − prod_lazy) mod m` in `[0, 2m)`:
    //   - `acc ≥ prod_lazy` ⇒ `diff = acc − prod_lazy ∈ [0, m)`
    //   - `acc <  prod_lazy` ⇒ `acc − prod_lazy + 2m ∈ [1, 2m − 1]`
    //
    // Barrett guarantees `2m < 2^BITS`, so the wrapping add of `2m` to the
    // wrap-borrowed `diff` produces the correct value modulo `2^BITS`.
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m = Simd::splat(modulus.value());
    let two_m = m + m;
    let (acc_chunks, acc_rem) = acc.as_chunks_mut::<N>();
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    for ((accc, ac), bc) in acc_chunks.iter_mut().zip(a_chunks).zip(b_chunks) {
        let accv = Simd::from_array(*accc);
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let prod_lazy = sm.lazy_reduce_mul(av, bv);
        let diff = accv - prod_lazy;
        *accc = accv
            .simd_lt(prod_lazy)
            .select(diff + two_m, diff)
            .to_array();
    }
    super::slice::scalar_lazy_reduce_sub_mul_slice_assign(modulus, acc_rem, a_rem, b_rem);
}

#[inline]
pub fn lazy_reduce_mul_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
    c: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let (a_chunks, a_rem) = a.as_chunks::<N>();
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (c_chunks, c_rem) = c.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for (((ac, bc), cc), oc) in a_chunks.iter().zip(b_chunks).zip(c_chunks).zip(o_chunks) {
        let av = Simd::from_array(*ac);
        let bv = Simd::from_array(*bc);
        let cv = Simd::from_array(*cc);
        *oc = sm.lazy_reduce_mul_add(av, bv, cv).to_array();
    }
    super::slice::scalar_lazy_reduce_mul_add_slice_to(modulus, a_rem, b_rem, c_rem, o_rem);
}

#[inline]
pub fn lazy_reduce_scalar_mul_add_slice_to<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    scalar: T,
    b: &[T],
    c: &[T],
    output: &mut [T],
) where
    Simd<T, N>: SimdArray<T, N>,
{
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let sv = Simd::splat(scalar);
    let (b_chunks, b_rem) = b.as_chunks::<N>();
    let (c_chunks, c_rem) = c.as_chunks::<N>();
    let (o_chunks, o_rem) = output.as_chunks_mut::<N>();
    for ((bc, cc), oc) in b_chunks.iter().zip(c_chunks).zip(o_chunks) {
        let bv = Simd::from_array(*bc);
        let cv = Simd::from_array(*cc);
        *oc = sm.lazy_reduce_mul_add(sv, bv, cv).to_array();
    }
    super::slice::scalar_lazy_reduce_scalar_mul_add_slice_to(modulus, scalar, b_rem, c_rem, o_rem);
}

// ---------------------------------------------------------------------------
// SIMD dot_product
//
// Outer chunk size = `K * N`, where `K = super::slice::DOT_PRODUCT_INNER_CHUNK`
// (currently 16). Inside each outer chunk we accumulate `K` SIMD widening
// products into a `[Simd<T, N>; 2]` double-word per lane, then collapse the
// double-word back into a single SIMD word in `[0, m)` via Barrett + the
// `min(v, v - m)` reduce_once trick. Cross-chunk accumulation stays in `[0, m)`
// lane-wise via `simd_reduce_add`, so the running SIMD accumulator never grows.
// After the chunked loop, a horizontal modular sum collapses the N lanes to a
// scalar, and any tail shorter than `K * N` is handled by the scalar helper.
//
// Hi-limb safety: each scalar widening product has `hi < m^2 / 2^BITS`, and the
// lo-limb's running sum can carry at most `K - 1` extra units into hi. With
// `m < 2^(BITS - 1)` enforced by `BarrettModulus::new` and `K ≤ 16`, both
// limbs stay strictly below `2^BITS` — identical bound to the scalar path.
// ---------------------------------------------------------------------------

#[inline]
pub fn reduce_dot_product<T: SimdUnsignedInteger, const N: usize>(
    modulus: BarrettModulus<T>,
    a: &[T],
    b: &[T],
) -> T
where
    Simd<T, N>: SimdArray<T, N>,
{
    assert_eq!(a.len(), b.len(), "reduce_dot_product: length mismatch");

    let k = super::slice::DOT_PRODUCT_INNER_CHUNK;
    let outer = k * N;

    let sm = SimdBarrettModulus::<T, N>::from(modulus);
    let m_simd = Simd::splat(modulus.value());
    let one = Simd::<T, N>::splat(T::ONE);
    let zero = Simd::<T, N>::splat(T::ZERO);

    let mut total_acc = Simd::<T, N>::splat(T::ZERO);

    let mut a_outer = a.chunks_exact(outer);
    let mut b_outer = b.chunks_exact(outer);

    for (a_chunk, b_chunk) in (&mut a_outer).zip(&mut b_outer) {
        // Each outer chunk is exactly `K * N` elements, so the inner
        // `as_chunks::<N>` always splits into `K` lane-wide subchunks with
        // an empty tail.
        let (a_lanes, _) = a_chunk.as_chunks::<N>();
        let (b_lanes, _) = b_chunk.as_chunks::<N>();
        let mut lo = Simd::<T, N>::splat(T::ZERO);
        let mut hi = Simd::<T, N>::splat(T::ZERO);
        for (a_n, b_n) in a_lanes.iter().zip(b_lanes) {
            let av = Simd::<T, N>::from_array(*a_n);
            let bv = Simd::<T, N>::from_array(*b_n);
            let (l, h) = av.widening_mul(bv);
            // Inline carry: lighter than `CarryingAdd::carrying_add` on
            // AVX2, where the latter expands to add + mask cast + sub +
            // 2x compare + or to surface a hi-carry we never use.
            let new_lo = lo + l;
            let carry = new_lo.simd_lt(lo).select(one, zero);
            lo = new_lo;
            hi = hi + h + carry;
        }
        let r = sm.lazy_reduce([lo, hi]);
        let r = simd_reduce_once(r, m_simd);
        total_acc = simd_reduce_add(total_acc, r, m_simd);
    }

    let tail_result =
        super::slice::scalar_reduce_dot_product(modulus, a_outer.remainder(), b_outer.remainder());

    let lanes = total_acc.to_array();
    let mut result = T::ZERO;
    for &v in lanes.iter() {
        result = modulus.reduce_add(result, v);
    }
    modulus.reduce_add(result, tail_result)
}
