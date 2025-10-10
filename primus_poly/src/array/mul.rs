use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::UintModulus;
use primus_reduce::ops::{ReduceAddAssign, ReduceMul, ReduceMulAdd, ReduceMulAssign};

use super::{ArrayBase, Data, DataMut, RawData};

impl<S, T> ArrayBase<S>
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
    pub fn add_mul_scalar_assign<M, A>(&mut self, rhs: &ArrayBase<A>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut()
            .zip(rhs)
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }

    #[inline]
    pub fn mul_element_wise_assign<M, A>(&mut self, rhs: &ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_mul_assign(a, b));
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
        rhs: &ArrayBase<A>,
        scalar: ShoupFactor<T>,
        modulus: T,
    ) where
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut().zip(rhs).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        });
    }
}

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn mul_element_wise_inplace<M, A>(&self, rhs: &Self, result: &mut ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_mul(a, b));
    }
}
