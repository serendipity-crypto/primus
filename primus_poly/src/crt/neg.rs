use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use crate::{ArrayBase, Data, DataMut, RawData};

use super::CrtPolynomial;

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg<M>(mut self, moduli: &[M], poly_length: usize) -> Self
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.neg_assign(moduli, poly_length);
        self
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, moduli: &[M], poly_length: usize)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.iter_mut(poly_length)
            .zip(moduli)
            .for_each(|(poly, modulus)| ArrayBase(poly).neg_assign(*modulus));
    }
}

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M, A>(
        &self,
        result: &mut CrtPolynomial<A, T>,
        moduli: &[M],
        poly_length: usize,
    ) where
        M: Copy + ReduceNeg<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(self.iter(poly_length), result.iter_mut(poly_length), moduli).for_each(
            |(xs, ys, modulus)| {
                ArrayBase(xs).neg_inplace(&mut ArrayBase(ys), *modulus);
            },
        );
    }
}
