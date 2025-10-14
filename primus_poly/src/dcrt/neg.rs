use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
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
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, moduli: &[M], poly_length: usize)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.iter_each_modulus_mut(poly_length)
            .zip(moduli)
            .for_each(|(poly, modulus)| ArrayBase(poly).neg_assign(*modulus));
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M, A>(
        &self,
        result: &mut DcrtPolynomial<A, T>,
        moduli: &[M],
        poly_length: usize,
    ) where
        M: Copy + ReduceNeg<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            result.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, modulus)| {
            ArrayBase(xs).neg_inplace(&mut ArrayBase(ys), *modulus);
        });
    }
}
