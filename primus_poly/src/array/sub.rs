use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceSub, ReduceSubAssign};

use super::{ArrayBase, Data, DataMut, DataOwned, RawData};

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub<M, A: RawData<Elem = T> + Data>(mut self, rhs: &ArrayBase<A>, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
    {
        self.sub_assign(rhs, modulus);
        self
    }
}

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_assign<M, A: RawData<Elem = T> + Data>(&mut self, rhs: &ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        izip!(self, rhs).for_each(|(a, &b)| modulus.reduce_sub_assign(a, b));
    }
}

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self - rhs` according to `modulus`.
    #[inline]
    pub fn sub_inplace<M, A: RawData<Elem = T> + DataMut>(
        &self,
        rhs: &Self,
        result: &mut ArrayBase<A>,
        modulus: M,
    ) where
        M: Copy + ReduceSub<T, Output = T>,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        debug_assert_eq!(self.0.len(), result.0.len());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_sub(a, b));
    }
}
