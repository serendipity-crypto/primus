use primus_integer::UnsignedInteger;
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use crate::{Data, DataMut, RawData};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<M, A>(mut self, rhs: &NttPolynomial<A, T>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M, A>(&mut self, rhs: &NttPolynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
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
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &NttPolynomial<A, T>,
        result: &mut NttPolynomial<B, T>,
        modulus: M,
    ) where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());
        self.iter()
            .zip(rhs.iter())
            .zip(result.iter_mut())
            .for_each(|((&a, &b), c)| *c = modulus.reduce_add(a, b))
    }
}
