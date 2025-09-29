use integer::{BigIntegerOps, UnsignedInteger, izip};

use super::BigUintPolynomial;

impl<T: UnsignedInteger> BigUintPolynomial<T> {
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub(mut self, rhs: &Self, modulus: &[T]) -> Self {
        self.sub_assign(rhs, modulus);
        self
    }

    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_assign(&mut self, rhs: &Self, modulus: &[T]) {
        debug_assert_eq!(self.big_uint_poly.len(), rhs.big_uint_poly.len());
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .zip(rhs.iter(value_len))
            .for_each(|(a, b)| {
                a.slice_sub_modulo_assign(b, modulus);
            });
    }

    /// Performs `result = self - rhs` according to `modulus`.
    #[inline]
    pub fn sub_inplace(&self, rhs: &Self, result: &mut Self, modulus: &[T]) {
        debug_assert_eq!(self.big_uint_poly.len(), rhs.big_uint_poly.len());
        debug_assert_eq!(self.big_uint_poly.len(), result.big_uint_poly.len());
        let value_len = modulus.len();
        izip!(
            self.iter(value_len),
            rhs.iter(value_len),
            result.iter_mut(value_len)
        )
        .for_each(|(a, b, c)| {
            a.slice_sub_modulo_inplace(b, c, modulus);
        });
    }
}
