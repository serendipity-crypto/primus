use primus_integer::{BigIntegerOps, UnsignedInteger, izip};

use crate::{Data, DataMut, RawData};

use super::BigUintPolynomial;

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<A>(mut self, rhs: &BigUintPolynomial<A>, modulus: &[T]) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<A>(&mut self, rhs: &BigUintPolynomial<A>, modulus: &[T])
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .zip(rhs.iter(value_len))
            .for_each(|(a, b)| {
                a.slice_add_modulo_assign(b, modulus);
            });
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<A, B>(
        &self,
        rhs: &BigUintPolynomial<A>,
        result: &mut BigUintPolynomial<B>,
        modulus: &[T],
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), rhs.len());
        debug_assert_eq!(self.len(), result.len());
        let value_len = modulus.len();
        izip!(
            self.iter(value_len),
            rhs.iter(value_len),
            result.iter_mut(value_len)
        )
        .for_each(|(a, b, c)| {
            a.slice_add_modulo_inplace(b, c, modulus);
        });
    }
}
