use primus_integer::{Data, DataMut, RawData, UnsignedInteger};
use primus_reduce::ops::{ReduceNeg, ReduceNegAssign};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
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

    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
    }
}

impl<S, T> NttPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs the unary `-` operation.
    #[inline]
    pub fn neg_inplace<M, A>(&self, result: &mut NttPolynomial<A>, modulus: M)
    where
        M: Copy + ReduceNeg<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.poly_length(), result.poly_length());
        self.iter()
            .zip(result.iter_mut())
            .for_each(|(&v, d)| *d = modulus.reduce_neg(v));
    }
}
