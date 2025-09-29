use primus_integer::{UnsignedInteger, izip};
use reduce::ops::{ReduceMul, ReduceMulAdd, ReduceMulAssign};

use crate::NttPolynomial;

use super::DcrtPolynomial;

impl<T: Copy> DcrtPolynomial<T> {
    /// Performs `self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, moduli: &[M]) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, moduli);
        self
    }

    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, moduli: &[M])
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_mut()
            .zip(moduli)
            .for_each(|(poly, modulus)| poly.mul_scalar_assign(scalar, *modulus))
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, moduli: &[M])
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        izip!(self, rhs, moduli).for_each(|(xs, ys, modulus)| {
            xs.add_mul_scalar_assign(ys, scalar, *modulus);
        });
    }
}

impl<T: UnsignedInteger> DcrtPolynomial<T> {
    /// Performs `self * rhs` according to `moduli`.
    #[inline]
    pub fn mul<M>(mut self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_assign(rhs, moduli);
        self
    }

    /// Performs `self *= rhs` according to `moduli`.
    #[inline]
    pub fn mul_assign<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: Copy + ReduceMulAssign<T>,
    {
        izip!(self, rhs, moduli).for_each(|(xs, ys, modulus)| xs.mul_assign(ys, *modulus))
    }

    /// Performs `result = self * rhs` according to `moduli`.
    #[inline]
    pub fn mul_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: Copy + ReduceMul<T, Output = T>,
    {
        izip!(self, rhs, result, moduli).for_each(
            |(xs, ys, zs, modulus): (
                &NttPolynomial<T>,
                &NttPolynomial<T>,
                &mut NttPolynomial<T>,
                &M,
            )| {
                xs.mul_inplace(ys, zs, *modulus);
            },
        )
    }
}
