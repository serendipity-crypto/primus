use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::UintModulus;
use primus_reduce::ops::*;

use crate::{Data, DataMut, DataOwned, RawData};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_assign(scalar, modulus);
        self
    }

    /// Performs `self * rhs` according to `modulus`.
    #[inline]
    pub fn mul<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_assign(rhs, modulus);
        self
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
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
    pub fn add_mul_scalar_assign<M, A>(&mut self, rhs: &NttPolynomial<A, T>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut()
            .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_factor_assign<A>(
        &mut self,
        rhs: &NttPolynomial<A, T>,
        scalar: ShoupFactor<T>,
        modulus: T,
    ) where
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut().zip(rhs.iter()).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        })
    }

    /// Performs `self *= rhs` according to `modulus`.
    #[inline]
    pub fn mul_assign<M, A>(&mut self, rhs: &NttPolynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(a, &b)| modulus.reduce_mul_assign(a, b));
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self * rhs` according to `modulus`.
    #[inline]
    pub fn mul_inplace<M, A>(&self, rhs: &Self, result: &mut NttPolynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(self.iter(), rhs.iter(), result.iter_mut())
            .for_each(|(&a, &b, c)| *c = modulus.reduce_mul(a, b));
    }
}
