use primus_integer::{BigIntegerOps, UnsignedInteger, izip};

use crate::{Data, DataMut, DataOwned, RawData};

use super::BigUintPolynomial;

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub<A>(mut self, rhs: &BigUintPolynomial<A>, modulus: &[T]) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign(rhs, modulus);
        self
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_assign<A>(&mut self, rhs: &BigUintPolynomial<A>, modulus: &[T])
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .zip(rhs.iter(value_len))
            .for_each(|(a, b)| {
                a.slice_sub_modulo_assign(b, modulus);
            });
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self - rhs` according to `modulus`.
    #[inline]
    pub fn sub_inplace<A>(&self, rhs: &Self, result: &mut BigUintPolynomial<A>, modulus: &[T])
    where
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        debug_assert_eq!(self.0.len(), result.0.len());
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
