use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign, ReduceMulAdd};

use super::{ArrayBase, Data, DataMut, RawData};

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `modulus`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &ArrayBase<A>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_element_wise_assign(rhs, modulus);
        self
    }

    /// Performs `self += rhs` according to `modulus`.
    #[inline]
    pub fn add_element_wise_assign<M, A>(&mut self, rhs: &ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
        izip!(self, rhs).for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }

    /// Performs `self += b * c` according to `modulus`.
    #[inline]
    pub fn add_mul_element_wise_assign<M, A, B>(
        &mut self,
        b: &ArrayBase<A>,
        c: &ArrayBase<B>,
        modulus: M,
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), b.len());
        izip!(self, b, c).for_each(|(a, &b, &c)| *a = modulus.reduce_mul_add(b, c, *a));
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `modulus`.
    #[inline]
    pub fn add_element_wise_inplace<M, A, B>(
        &self,
        rhs: &ArrayBase<A>,
        result: &mut ArrayBase<B>,
        modulus: M,
    ) where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), rhs.len());
        debug_assert_eq!(self.len(), result.len());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_add(a, b));
    }
}
