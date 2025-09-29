use std::simd::{LaneCount, Simd, SupportedLaneCount};

use primus_integer::{CarryingAdd, CarryingMul, SimdArray, SimdUnsignedInteger, WideningMul};
use reduce::lazy_ops::*;

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
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    value: Simd<T, N>,
    ratio: [Simd<T, N>; 2],
}

impl<T: SimdUnsignedInteger, const N: usize> From<BarrettModulus<T>> for SimdBarrettModulus<T, N>
where
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_reduce_assign(self, value: &mut Simd<T, N>) {
        *value = self.lazy_reduce(*value);
    }
}

impl<T: SimdUnsignedInteger, const N: usize> LazyReduceMul<Simd<T, N>> for SimdBarrettModulus<T, N>
where
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
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
    LaneCount<N>: SupportedLaneCount,
    Simd<T, N>: SimdArray<T, N>,
{
    #[inline]
    fn lazy_reduce_mul_add_assign(self, a: &mut Simd<T, N>, b: Simd<T, N>, c: Simd<T, N>) {
        *a = self.lazy_reduce(a.carrying_mul(b, c));
    }
}
