use primus_integer::izip;
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use super::CrtPolynomial;

impl<T: Copy> CrtPolynomial<T> {
    /// Performs `self + rhs` according to `moduli`.
    #[inline]
    pub fn add<M>(mut self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_assign(rhs, moduli);
        self
    }

    /// Performs `self += rhs` according to `moduli`.
    #[inline]
    pub fn add_assign<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: Copy + ReduceAddAssign<T>,
    {
        izip!(self, rhs, moduli).for_each(|(xs, ys, modulus)| {
            xs.add_assign(ys, *modulus);
        });
    }

    /// Performs `result = self + rhs` according to `moduli`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        izip!(self, rhs, result, moduli).for_each(|(xs, ys, zs, modulus)| {
            xs.add_inplace(ys, zs, *modulus);
        });
    }
}
