use primus_integer::{Data, DataMut, RawData, UnsignedInteger, izip};
use primus_reduce::ops::{ReduceSub, ReduceSubAssign};

use super::ArrayBase;

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &ArrayBase<A>, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_element_wise_assign(rhs, modulus);
        self
    }

    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_element_wise_assign<M, A>(&mut self, rhs: &ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.len(), rhs.len());
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
    pub fn sub_element_wise_inplace<M, A, B>(
        &self,
        rhs: &ArrayBase<A>,
        result: &mut ArrayBase<B>,
        modulus: M,
    ) where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), rhs.len());
        debug_assert_eq!(self.len(), result.len());
        izip!(self, rhs, result).for_each(|(&a, &b, c)| *c = modulus.reduce_sub(a, b));
    }

    /// Performs `rhs = self - rhs` according to `modulus`.
    #[inline]
    pub fn sub_element_wise_to_right<M, A>(&self, rhs: &mut ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.len(), rhs.len());
        izip!(self, rhs).for_each(|(&a, b)| *b = modulus.reduce_sub(a, *b));
    }
}
