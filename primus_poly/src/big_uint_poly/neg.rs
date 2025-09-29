use integer::{BigIntegerOps, UnsignedInteger};

use super::BigUintPolynomial;

impl<T: UnsignedInteger> BigUintPolynomial<T> {
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg(mut self, modulus: &[T]) -> Self {
        self.neg_assign(modulus);
        self
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign(&mut self, modulus: &[T]) {
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .for_each(|v| v.slice_neg_modulo_assign(modulus));
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace(&self, result: &mut Self, modulus: &[T]) {
        debug_assert_eq!(self.big_uint_poly.len(), result.big_uint_poly.len());
        let value_len = modulus.len();
        result
            .iter_mut(value_len)
            .zip(self.iter(value_len))
            .for_each(|(d, v)| v.slice_neg_modulo_inplace(d, modulus));
    }
}
