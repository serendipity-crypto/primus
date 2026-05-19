//! Scalar trait impls — all delegate to `Barrett50Modulus::inner`
//! (`BarrettModulus<u64>`). Scalar paths are not IFMA-accelerated; the
//! win lives in the SIMD slice kernels in [`super::simd_ifma`].

use primus_integer::UnsignedInteger;
use primus_reduce::{ReduceError, lazy_ops::*, ops::*};

use super::Barrett50Modulus;

// ---------------------------------------------------------------------------
// LazyReduce family
// ---------------------------------------------------------------------------

impl LazyReduce<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce(self, value: u64) -> u64 {
        self.inner.lazy_reduce(value)
    }
}

impl LazyReduce<[u64; 2]> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce(self, value: [u64; 2]) -> u64 {
        self.inner.lazy_reduce(value)
    }
}

impl LazyReduce<(u64, u64)> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce(self, value: (u64, u64)) -> u64 {
        self.inner.lazy_reduce(value)
    }
}

impl LazyReduce<&[u64]> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce(self, value: &[u64]) -> u64 {
        self.inner.lazy_reduce(value)
    }
}

impl LazyReduceAssign<u64> for Barrett50Modulus {
    #[inline]
    fn lazy_reduce_assign(self, value: &mut u64) {
        self.inner.lazy_reduce_assign(value);
    }
}

impl LazyReduceMul<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce_mul(self, a: u64, b: u64) -> u64 {
        self.inner.lazy_reduce_mul(a, b)
    }
}

impl LazyReduceMulAssign<u64> for Barrett50Modulus {
    #[inline]
    fn lazy_reduce_mul_assign(self, a: &mut u64, b: u64) {
        self.inner.lazy_reduce_mul_assign(a, b);
    }
}

impl LazyReduceMulAdd<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn lazy_reduce_mul_add(self, a: u64, b: u64, c: u64) -> u64 {
        self.inner.lazy_reduce_mul_add(a, b, c)
    }
}

impl LazyReduceMulAddAssign<u64> for Barrett50Modulus {
    #[inline]
    fn lazy_reduce_mul_add_assign(self, a: &mut u64, b: u64, c: u64) {
        self.inner.lazy_reduce_mul_add_assign(a, b, c);
    }
}

// ---------------------------------------------------------------------------
// Reduce family
// ---------------------------------------------------------------------------

impl Reduce<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce(self, value: u64) -> u64 {
        self.inner.reduce(value)
    }
}

impl Reduce<[u64; 2]> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce(self, value: [u64; 2]) -> u64 {
        self.inner.reduce(value)
    }
}

impl Reduce<(u64, u64)> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce(self, value: (u64, u64)) -> u64 {
        self.inner.reduce(value)
    }
}

impl Reduce<&[u64]> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce(self, value: &[u64]) -> u64 {
        self.inner.reduce(value)
    }
}

impl ReduceAssign<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_assign(self, value: &mut u64) {
        self.inner.reduce_assign(value);
    }
}

impl ReduceOnce<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_once(self, value: u64) -> u64 {
        self.inner.reduce_once(value)
    }
}

impl ReduceOnceAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_once_assign(self, value: &mut u64) {
        self.inner.reduce_once_assign(value);
    }
}

impl ReduceAdd<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_add(self, a: u64, b: u64) -> u64 {
        self.inner.reduce_add(a, b)
    }
}

impl ReduceAddAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_add_assign(self, a: &mut u64, b: u64) {
        self.inner.reduce_add_assign(a, b);
    }
}

impl ReduceDouble<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_double(self, value: u64) -> u64 {
        self.inner.reduce_double(value)
    }
}

impl ReduceDoubleAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_double_assign(self, value: &mut u64) {
        self.inner.reduce_double_assign(value);
    }
}

impl ReduceSub<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_sub(self, a: u64, b: u64) -> u64 {
        self.inner.reduce_sub(a, b)
    }
}

impl ReduceSubAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_sub_assign(self, a: &mut u64, b: u64) {
        self.inner.reduce_sub_assign(a, b);
    }
}

impl ReduceNeg<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_neg(self, value: u64) -> u64 {
        self.inner.reduce_neg(value)
    }
}

impl ReduceNegAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_neg_assign(self, value: &mut u64) {
        self.inner.reduce_neg_assign(value);
    }
}

impl ReduceMul<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn reduce_mul(self, a: u64, b: u64) -> u64 {
        self.inner.reduce_mul(a, b)
    }
}

impl ReduceMulAssign<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_mul_assign(self, a: &mut u64, b: u64) {
        self.inner.reduce_mul_assign(a, b);
    }
}

impl ReduceSquare<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn reduce_square(self, value: u64) -> u64 {
        self.inner.reduce_square(value)
    }
}

impl ReduceSquareAssign<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_square_assign(self, value: &mut u64) {
        self.inner.reduce_square_assign(value);
    }
}

impl ReduceMulAdd<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn reduce_mul_add(self, a: u64, b: u64, c: u64) -> u64 {
        self.inner.reduce_mul_add(a, b, c)
    }
}

impl ReduceMulAddAssign<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_mul_add_assign(self, a: &mut u64, b: u64, c: u64) {
        self.inner.reduce_mul_add_assign(a, b, c);
    }
}

impl TryReduceInv<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn try_reduce_inv(self, value: u64) -> Result<u64, ReduceError<u64>> {
        self.inner.try_reduce_inv(value)
    }
}

impl ReduceInv<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline(always)]
    fn reduce_inv(self, value: u64) -> u64 {
        self.inner.reduce_inv(value)
    }
}

impl ReduceInvAssign<u64> for Barrett50Modulus {
    #[inline(always)]
    fn reduce_inv_assign(self, value: &mut u64) {
        self.inner.reduce_inv_assign(value);
    }
}

impl ReduceDiv<u64> for Barrett50Modulus {
    type Output = u64;
    #[inline]
    fn reduce_div(self, a: u64, b: u64) -> u64 {
        self.inner.reduce_div(a, b)
    }
}

impl ReduceDivAssign<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_div_assign(self, a: &mut u64, b: u64) {
        self.inner.reduce_div_assign(a, b);
    }
}

impl ReduceExp<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_exp<E: UnsignedInteger>(self, base: u64, exp: E) -> u64 {
        self.inner.reduce_exp(base, exp)
    }
}

impl ReduceExpPowOf2<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_exp_power_of_2(self, base: u64, exp_log: u32) -> u64 {
        self.inner.reduce_exp_power_of_2(base, exp_log)
    }
}
