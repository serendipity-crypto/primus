use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::UintModulus;
use primus_reduce::ops::{ReduceAddAssign, ReduceMul, ReduceMulAdd, ReduceMulAssign};

use super::{Array, ArrayMut, ArrayRef};

impl<T: Copy> Array<T> {}

impl<'a, T: Copy> ArrayMut<'a, T> {
    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.0
            .iter_mut()
            .for_each(|v| modulus.reduce_mul_assign(v, scalar))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: ArrayRef<'_, T>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.0
            .iter_mut()
            .zip(rhs.0)
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }

    #[inline]
    pub fn mul_element_wise_assign<M>(&mut self, rhs: ArrayRef<'_, T>, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.0
            .iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_mul_assign(a, b));
    }
}

impl<'a, T: UnsignedInteger> ArrayMut<'a, T> {
    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.0
            .iter_mut()
            .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_factor_assign(
        &mut self,
        rhs: ArrayRef<'_, T>,
        scalar: ShoupFactor<T>,
        modulus: T,
    ) {
        self.0.iter_mut().zip(rhs).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        });
    }
}

impl<'a, T: Copy> ArrayRef<'a, T> {
    #[inline]
    pub fn mul_element_wise_inplace<M>(
        self,
        rhs: ArrayRef<'_, T>,
        result: &mut ArrayMut<'_, T>,
        modulus: M,
    ) where
        M: Copy + ReduceMul<T, Output = T>,
    {
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_mul(a, b));
    }
}
