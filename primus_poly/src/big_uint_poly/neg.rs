use primus_integer::{BigUint, Data, DataMut, RawData, UnsignedInteger};

use super::BigUintPolynomial;

impl<S, T> BigUintPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg<A>(mut self, modulus: &BigUint<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        self.neg_assign(modulus);
        self
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<A>(&mut self, modulus: &BigUint<A>)
    where
        A: RawData<Elem = T> + Data,
    {
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .for_each(|mut v| v.neg_modulo_assign(modulus));
    }
}

impl<S, T> BigUintPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<A, B>(&self, result: &mut BigUintPolynomial<A>, modulus: &BigUint<B>)
    where
        A: RawData<Elem = T> + DataMut,
        B: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), result.len());
        let value_len = modulus.len();
        result
            .iter_mut(value_len)
            .zip(self.iter(value_len))
            .for_each(|(mut d, v)| v.neg_modulo_inplace(&mut d, modulus));
    }
}
