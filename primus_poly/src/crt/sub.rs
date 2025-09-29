use primus_integer::izip;
use reduce::ops::{ReduceSub, ReduceSubAssign};

use super::CrtPolynomial;

impl<T: Copy> CrtPolynomial<T> {
    /// Performs `self - rhs` according to `moduli`.
    #[inline]
    pub fn sub<M>(mut self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
    {
        self.sub_assign(rhs, moduli);
        self
    }

    /// Performs `self -= rhs` according to `moduli`.
    #[inline]
    pub fn sub_assign<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: Copy + ReduceSubAssign<T>,
    {
        izip!(self, rhs, moduli).for_each(|(xs, ys, modulus)| {
            xs.sub_assign(ys, *modulus);
        });
    }

    /// Performs `result = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: Copy + ReduceSub<T, Output = T>,
    {
        izip!(self, rhs, result, moduli).for_each(|(xs, ys, zs, modulus)| {
            xs.sub_inplace(ys, zs, *modulus);
        });
    }
}
