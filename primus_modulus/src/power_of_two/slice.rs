//! Scalar slice impls for [`PowOf2Modulus`].
//!
//! Each operation reduces by masking with `self.mask`, which is just bit-AND
//! and lets the compiler auto-vectorize most loops on its own. The impls
//! below loop over the slice using the existing scalar `Reduce*` impls.

use primus_integer::UnsignedInteger;
use primus_reduce::{ReduceMulAdd, ReduceMulAddSlice};

use super::PowOf2Modulus;

impl<T: UnsignedInteger> ReduceMulAddSlice<T> for PowOf2Modulus<T> {
    #[inline]
    fn reduce_add_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]) {
        debug_assert_eq!(acc.len(), a.len());
        debug_assert_eq!(acc.len(), b.len());
        acc.iter_mut()
            .zip(a)
            .zip(b)
            .for_each(|((acc, &a), &b)| *acc = self.reduce_mul_add(a, b, *acc));
    }

    #[inline]
    fn reduce_sub_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]) {
        debug_assert_eq!(acc.len(), a.len());
        debug_assert_eq!(acc.len(), b.len());
        let mask = self.mask();
        acc.iter_mut()
            .zip(a)
            .zip(b)
            .for_each(|((acc, &a), &b)| *acc = (*acc).wrapping_sub(a.wrapping_mul(b)) & mask);
    }

    #[inline]
    fn reduce_mul_add_slice_to(self, a: &[T], b: &[T], c: &[T], output: &mut [T]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), c.len());
        debug_assert_eq!(a.len(), output.len());
        a.iter()
            .zip(b)
            .zip(c)
            .zip(output)
            .for_each(|(((&a, &b), &c), o)| *o = self.reduce_mul_add(a, b, c));
    }

    #[inline]
    fn reduce_scalar_mul_add_slice_to(self, scalar: T, b: &[T], c: &[T], output: &mut [T]) {
        debug_assert_eq!(b.len(), c.len());
        debug_assert_eq!(b.len(), output.len());
        b.iter()
            .zip(c)
            .zip(output)
            .for_each(|((&b, &c), o)| *o = self.reduce_mul_add(scalar, b, c));
    }
}
