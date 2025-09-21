use integer::UnsignedInteger;
use num_traits::{ConstZero, One, Zero};
use primus_factor::{FactorMul, ShoupFactor};
use reduce::ops::*;
use uint_modulus::UintModulus;

use super::NttPolynomial;

impl<T> NttPolynomial<T>
where
    T: Copy + Zero + ConstZero + One + PartialEq,
{
    /// Multiply `self` with a scalar.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Multiply `self` with a scalar assign.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        if scalar.is_zero() {
            self.set_zero();
        } else if !scalar.is_one() {
            self.iter_mut()
                .for_each(|v| modulus.reduce_mul_assign(v, scalar))
        }
    }

    /// Add the multiply result `rhs` with a scalar inplace.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, modulus: M)
    where
        M: Copy + ReduceAddAssign<T> + ReduceMulAdd<T, Output = T>,
    {
        if scalar.is_one() {
            self.add_assign(rhs, modulus);
        } else if !scalar.is_zero() {
            self.iter_mut()
                .zip(rhs.iter())
                .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
        }
    }
}

impl<T: UnsignedInteger> NttPolynomial<T> {
    /// Multiply `self` with a shoup scalar.
    #[inline]
    pub fn mul_factor_scalar(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_scalar_assign(scalar, modulus);
        self
    }

    /// Multiply `self` with a shoup scalar assign.
    #[inline]
    pub fn mul_factor_scalar_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        if scalar.value().is_zero() {
            self.set_zero();
        } else if !scalar.value().is_one() {
            self.iter_mut()
                .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
        }
    }

    /// Multiply `self` with a scalar and add to self.
    #[inline]
    pub fn add_mul_factor_scalar_assign(&mut self, rhs: &Self, scalar: ShoupFactor<T>, modulus: T) {
        if scalar.value().is_one() {
            self.add_assign(rhs, UintModulus(modulus));
        } else if !scalar.value().is_zero() {
            self.iter_mut().zip(rhs.iter()).for_each(|(r, &v)| {
                UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
            })
        }
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

    /// Performs `destination = self * rhs` according to `modulus`.
    #[inline]
    pub fn mul_inplace<M>(&self, rhs: &Self, modulus: M, destination: &mut Self)
    where
        M: Copy + ReduceMul<T, Output = T>,
    {
        self.iter()
            .zip(rhs)
            .zip(destination)
            .for_each(|((&a, &b), c)| *c = modulus.reduce_mul(a, b));
    }
}
