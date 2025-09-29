use primus_integer::izip;
use reduce::ops::{ReduceSub, ReduceSubAssign};

use super::NttPolynomial;

impl<T: Copy> NttPolynomial<T> {
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
    {
        self.sub_assign(rhs, modulus);
        self
    }

    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_assign<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        self.iter_mut()
            .zip(rhs)
            .for_each(|(a, &b)| modulus.reduce_sub_assign(a, b));
    }

    /// Performs `result = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: Copy + ReduceSub<T, Output = T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());

        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_sub(a, b))
    }
}
