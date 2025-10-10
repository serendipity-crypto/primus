use primus_integer::UnsignedInteger;
use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use super::{ArrayBase, Data, DataMut, DataOwned, RawData};

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg<M>(mut self, modulus: M) -> Self
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.neg_assign(modulus);
        self
    }
}

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
    }
}

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M, A>(&self, result: &mut ArrayBase<A>, modulus: M)
    where
        M: Copy + ReduceNeg<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.0.len(), result.0.len());
        self.iter()
            .zip(result)
            .for_each(|(&v, d)| *d = modulus.reduce_neg(v));
    }
}
