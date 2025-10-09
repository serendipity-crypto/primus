use primus_integer::izip;
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use super::{Array, ArrayMut, ArrayRef};

impl<T: Copy> Array<T> {
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<M>(mut self, rhs: ArrayRef<'_, T>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M>(&mut self, rhs: ArrayRef<'_, T>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.to_mut().add_assign(rhs, modulus);
    }

    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        self.to_ref()
            .add_inplace(rhs.to_ref(), &mut result.to_mut(), modulus);
    }
}

impl<'a, T: Copy> ArrayMut<'a, T> {
    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M>(&mut self, rhs: ArrayRef<'_, T>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        izip!(self, rhs).for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }
}

impl<'a, T: Copy> ArrayRef<'a, T> {
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<M>(self, rhs: ArrayRef<'_, T>, result: &mut ArrayMut<'_, T>, modulus: M)
    where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        debug_assert_eq!(self.0.len(), result.0.len());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_add(a, b));
    }
}
