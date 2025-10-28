use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::UintModulus;
use primus_reduce::ops::*;

use super::{ArrayBase, Data, DataMut, RawData};

impl<S, T> ArrayBase<S, T>
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
            .for_each(|a| modulus.reduce_mul_assign(a, scalar))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_scalar_assign<M, A>(&mut self, rhs: &ArrayBase<A>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| *a = modulus.reduce_mul_add(b, scalar, *a));
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut()
            .for_each(|a| *a = scalar.factor_mul_modulo(*a, modulus))
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
        debug_assert_eq!(self.len(), rhs.len());
        self.iter_mut().zip(rhs).for_each(|(a, &b)| {
            UintModulus(modulus).reduce_add_assign(a, scalar.factor_mul_modulo(b, modulus))
        });
    }

    #[inline]
    pub fn mul_element_wise_assign<M, A>(&mut self, rhs: &ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_mul_assign(a, b));
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs element wise modular multiplication operation `result = self * rhs` according to `modulus`.
    #[inline]
    pub fn mul_element_wise_inplace<M, A, B>(
        &self,
        rhs: &ArrayBase<A>,
        result: &mut ArrayBase<B>,
        modulus: M,
    ) where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), rhs.len());
        debug_assert_eq!(self.len(), result.len());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_mul(a, b));
    }

    /// Performs `result = scalar * self` according to `modulus`.
    #[inline]
    pub fn mul_scalar_inplace<M, A>(&self, scalar: T, result: &mut ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), result.len());
        self.iter()
            .zip(result.iter_mut())
            .for_each(|(&a, b)| *b = modulus.reduce_mul(a, scalar));
    }

    /// Performs `result = scalar * self` according to `modulus`.
    #[inline]
    pub fn mul_factor_inplace<A>(
        &self,
        scalar: ShoupFactor<T>,
        result: &mut ArrayBase<A>,
        modulus: T,
    ) where
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), result.len());
        self.iter()
            .zip(result.iter_mut())
            .for_each(|(&a, b)| *b = scalar.factor_mul_modulo(a, modulus));
    }
}
