use integer::izip;
use reduce::ops::{ReduceNeg, ReduceNegAssign};

use super::CrtPolynomial;

impl<T: Copy> CrtPolynomial<T> {
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg<M>(mut self, moduli: &[M]) -> Self
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.neg_assign(moduli);
        self
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, moduli: &[M])
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.iter_mut()
            .zip(moduli)
            .for_each(|(poly, modulus)| poly.neg_assign(*modulus));
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M>(&self, result: &mut Self, moduli: &[M])
    where
        M: Copy + ReduceNeg<T, Output = T>,
    {
        izip!(self, result, moduli).for_each(|(xs, ys, modulus)| {
            xs.neg_inplace(ys, *modulus);
        });
    }
}
