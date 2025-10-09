use primus_integer::izip;
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use super::Polynomial;

impl<T: Copy> Polynomial<T> {
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }

    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_add(a, b));
    }
}
