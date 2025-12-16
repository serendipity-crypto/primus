use primus_integer::{BigUint, Data, DataMut, RawData, UnsignedInteger, izip};

use super::BigUintPolynomial;

impl<S, T> BigUintPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<A, B>(mut self, rhs: &BigUintPolynomial<A>, modulus: &BigUint<B>) -> Self
    where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<A, B>(&mut self, rhs: &BigUintPolynomial<A>, modulus: &BigUint<B>)
    where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        let value_len = modulus.len();
        self.iter_mut(value_len)
            .zip(rhs.iter(value_len))
            .for_each(|(mut a, b)| {
                a.add_modulo_assign(&b, modulus);
            });
    }
}

impl<S, T> BigUintPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<A, B, C>(
        &self,
        rhs: &BigUintPolynomial<A>,
        result: &mut BigUintPolynomial<B>,
        modulus: &BigUint<C>,
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
        C: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        debug_assert_eq!(self.len(), result.len());
        let value_len = modulus.len();
        izip!(
            self.iter(value_len),
            rhs.iter(value_len),
            result.iter_mut(value_len)
        )
        .for_each(|(a, b, mut c)| {
            a.add_modulo_inplace(&b, &mut c, modulus);
        });
    }
}
