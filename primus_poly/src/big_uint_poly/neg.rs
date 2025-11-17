use primus_integer::{BigIntegerOps, UnsignedInteger};

use crate::{Data, DataMut, RawData};

use super::BigUintPolynomial;

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
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
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<A>(&self, result: &mut BigUintPolynomial<A>, modulus: &[T])
    where
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), result.len());
        let value_len = modulus.len();
        result
            .iter_mut(value_len)
            .zip(self.iter(value_len))
            .for_each(|(d, v)| v.slice_neg_modulo_inplace(d, modulus));
    }
}
