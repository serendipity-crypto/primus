use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use crate::{Data, DataMut, RawData};

use super::Polynomial;

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add<M, A: RawData<Elem = T> + Data>(mut self, rhs: &Polynomial<A, T>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_assign<M, A>(&mut self, rhs: &Polynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut Polynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());
        izip!(&self.0, &rhs.0, &mut result.0).for_each(|(&a, &b, c)| *c = modulus.reduce_add(a, b));
    }
}
