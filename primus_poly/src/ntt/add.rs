use primus_integer::UnsignedInteger;
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use crate::{Data, DataMut, DataOwned, RawData};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<M, A: RawData<Elem = T> + Data>(
        mut self,
        rhs: &NttPolynomial<A, T>,
        modulus: M,
    ) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_assign(rhs, modulus);
        self
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M, A: RawData<Elem = T> + Data>(
        &mut self,
        rhs: &NttPolynomial<A, T>,
        modulus: M,
    ) where
        M: Copy + ReduceAddAssign<T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<M, A: RawData<Elem = T> + DataMut>(
        &self,
        rhs: &Self,
        result: &mut NttPolynomial<A, T>,
        modulus: M,
    ) where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());
        self.iter()
            .zip(rhs.iter())
            .zip(result.iter_mut())
            .for_each(|((&a, &b), c)| *c = modulus.reduce_add(a, b))
    }
}
