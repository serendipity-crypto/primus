use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use super::{Array, ArrayMut, ArrayRef};

impl<T: Copy> Array<T> {
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg<M>(mut self, modulus: M) -> Self
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.neg_assign(modulus);
        self
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.to_mut().neg_assign(modulus);
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M>(&self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceNeg<T, Output = T>,
    {
        self.to_ref().neg_inplace(&mut result.to_mut(), modulus);
    }
}

impl<'a, T: Copy> ArrayMut<'a, T> {
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.0.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
    }
}

impl<'a, T: Copy> ArrayRef<'a, T> {
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M>(self, result: &mut ArrayMut<'_, T>, modulus: M)
    where
        M: Copy + ReduceNeg<T, Output = T>,
    {
        debug_assert_eq!(self.0.len(), result.0.len());
        result
            .0
            .iter_mut()
            .zip(self)
            .for_each(|(d, &v)| *d = modulus.reduce_neg(v));
    }
}
