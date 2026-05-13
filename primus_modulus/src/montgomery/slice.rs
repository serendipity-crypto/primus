//! Scalar slice impls for [`MontgomeryModulus`].
//!
//! These loop over slices using the existing scalar `Reduce*` /
//! `LazyReduce*` impls on [`MontgomeryModulus`]. No SIMD path yet.

use primus_integer::UnsignedInteger;
use primus_reduce::{
    LazyReduceMulAdd, LazyReduceMulAddSlice, ReduceMulAdd, ReduceMulAddSlice, ReduceSubAssign,
};

use super::MontgomeryModulus;

impl<T: UnsignedInteger> ReduceMulAddSlice<T> for MontgomeryModulus<T> {
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
        // No scalar `reduce_sub_mul` exists; compute `a*b mod m` then sub.
        // `reduce_mul_add(a, b, 0)` gives the canonical product.
        acc.iter_mut().zip(a).zip(b).for_each(|((acc, &a), &b)| {
            let prod = self.reduce_mul_add(a, b, T::ZERO);
            self.reduce_sub_assign(acc, prod);
        });
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
    fn reduce_scalar_mul_add_slice_to(
        self,
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
            .for_each(|((&b, &c), o)| *o = self.reduce_mul_add(scalar, b, c));
    }
}

impl<T: UnsignedInteger> LazyReduceMulAddSlice<T> for MontgomeryModulus<T> {
    #[inline]
    fn lazy_reduce_add_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]) {
        debug_assert_eq!(acc.len(), a.len());
        debug_assert_eq!(acc.len(), b.len());
        acc.iter_mut()
            .zip(a)
            .zip(b)
            .for_each(|((acc, &a), &b)| *acc = self.lazy_reduce_mul_add(a, b, *acc));
    }

    #[inline]
    fn lazy_reduce_sub_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]) {
        debug_assert_eq!(acc.len(), a.len());
        debug_assert_eq!(acc.len(), b.len());
        acc.iter_mut().zip(a).zip(b).for_each(|((acc, &a), &b)| {
            let prod = self.lazy_reduce_mul_add(a, b, T::ZERO);
            self.reduce_sub_assign(acc, prod);
        });
    }

    #[inline]
    fn lazy_reduce_mul_add_slice_to(self, a: &[T], b: &[T], c: &[T], output: &mut [T]) {
        debug_assert_eq!(a.len(), b.len());
        debug_assert_eq!(a.len(), c.len());
        debug_assert_eq!(a.len(), output.len());
        a.iter()
            .zip(b)
            .zip(c)
            .zip(output)
            .for_each(|(((&a, &b), &c), o)| *o = self.lazy_reduce_mul_add(a, b, c));
    }

    #[inline]
    fn lazy_reduce_scalar_mul_add_slice_to(
        self,
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
            .for_each(|((&b, &c), o)| *o = self.lazy_reduce_mul_add(scalar, b, c));
    }
}
