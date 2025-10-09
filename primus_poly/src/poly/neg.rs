use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use super::Polynomial;

impl<T: Copy> Polynomial<T> {
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
        self.data
            .iter_mut()
            .for_each(|v| modulus.reduce_neg_assign(v));
    }

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M>(&self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceNeg<T, Output = T>,
    {
        debug_assert_eq!(self.poly_length(), result.poly_length());
        result
            .iter_mut()
            .zip(self)
            .for_each(|(d, &v)| *d = modulus.reduce_neg(v));
    }
}
