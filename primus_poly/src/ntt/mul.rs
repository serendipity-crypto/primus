use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::UintModulus;
use primus_reduce::ops::*;

use super::NttPolynomial;

impl<T: Copy> NttPolynomial<T> {
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_mut()
            .for_each(|v| modulus.reduce_mul_assign(v, scalar))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.iter_mut()
            .zip(rhs)
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }
}

impl<T: UnsignedInteger> NttPolynomial<T> {
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_assign(scalar, modulus);
        self
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut()
            .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_factor_assign(&mut self, rhs: &Self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut().zip(rhs.iter()).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        })
    }
}

impl<T: UnsignedInteger> NttPolynomial<T> {
    /// Performs `self * rhs` according to `modulus`.
    #[inline]
    pub fn mul<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_assign(rhs, modulus);
        self
    }

    /// Performs `self *= rhs` according to `modulus`.
    #[inline]
    pub fn mul_assign<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_mul_assign(a, b));
    }

    /// Performs `result = self * rhs` according to `modulus`.
    #[inline]
    pub fn mul_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceMul<T, Output = T>,
    {
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_mul(a, b));
    }
}
