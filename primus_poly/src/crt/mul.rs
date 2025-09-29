use primus_integer::izip;
use primus_reduce::ops::{ReduceMulAdd, ReduceMulAssign};

use super::CrtPolynomial;

impl<T: Copy> CrtPolynomial<T> {
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
