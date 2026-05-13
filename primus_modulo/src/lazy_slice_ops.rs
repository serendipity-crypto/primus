//! Value-side mirror of [`primus_reduce::lazy_slice_ops`].
//!
//! Results are in `[0, 2 * modulus)`; callers must perform a final
//! reduction (e.g. via [`crate::slice_ops::OnceModuloSlice`]) when a
//! canonical representative is required.

use primus_reduce::lazy_slice_ops::*;

/// Value-side mirror of [`LazyReduceMulSlice`].
pub trait LazyMulModuloSlice<M, B: ?Sized = Self> {
    /// `self[i] = self[i] * b[i] (mod 2 * modulus)` element-wise.
    fn lazy_mul_modulo_slice_assign(&mut self, b: &B, modulus: M);

    /// `output[i] = self[i] * b[i] (mod 2 * modulus)`.
    fn lazy_mul_modulo_slice_to(&self, b: &B, output: &mut Self, modulus: M);
}

impl<T, M, B> LazyMulModuloSlice<M, [B]> for [T]
where
    M: LazyReduceMulSlice<T, B> + Copy,
{
    #[inline(always)]
    fn lazy_mul_modulo_slice_assign(&mut self, b: &[B], modulus: M) {
        modulus.lazy_reduce_mul_slice_assign(self, b);
    }

    #[inline(always)]
    fn lazy_mul_modulo_slice_to(&self, b: &[B], output: &mut [T], modulus: M) {
        modulus.lazy_reduce_mul_slice_to(self, b, output);
    }
}

/// Value-side mirror of [`LazyReduceMulAddSlice`].
pub trait LazyMulAddModuloSlice<M, T> {
    /// `self[i] += a[i] * b[i] (mod 2 * modulus)`.
    fn lazy_add_mul_modulo_slice_assign(&mut self, a: &[T], b: &[T], modulus: M);

    /// `self[i] -= a[i] * b[i] (mod 2 * modulus)`.
    fn lazy_sub_mul_modulo_slice_assign(&mut self, a: &[T], b: &[T], modulus: M);

    /// `output[i] = self[i] * b[i] + c[i] (mod 2 * modulus)`.
    fn lazy_mul_add_modulo_slice_to(&self, b: &[T], c: &[T], output: &mut [T], modulus: M);

    /// `output[i] = scalar * self[i] + c[i] (mod 2 * modulus)`.
    fn lazy_scalar_mul_add_modulo_slice_to(&self, scalar: T, c: &[T], output: &mut [T], modulus: M);
}

impl<T, M> LazyMulAddModuloSlice<M, T> for [T]
where
    M: LazyReduceMulAddSlice<T> + Copy,
{
    #[inline(always)]
    fn lazy_add_mul_modulo_slice_assign(&mut self, a: &[T], b: &[T], modulus: M) {
        modulus.lazy_reduce_add_mul_slice_assign(self, a, b);
    }

    #[inline(always)]
    fn lazy_sub_mul_modulo_slice_assign(&mut self, a: &[T], b: &[T], modulus: M) {
        modulus.lazy_reduce_sub_mul_slice_assign(self, a, b);
    }

    #[inline(always)]
    fn lazy_mul_add_modulo_slice_to(&self, b: &[T], c: &[T], output: &mut [T], modulus: M) {
        modulus.lazy_reduce_mul_add_slice_to(self, b, c, output);
    }

    #[inline(always)]
    fn lazy_scalar_mul_add_modulo_slice_to(
        &self,
        scalar: T,
        c: &[T],
        output: &mut [T],
        modulus: M,
    ) {
        modulus.lazy_reduce_scalar_mul_add_slice_to(scalar, self, c, output);
    }
}
