use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceSub, ReduceSubAssign};

use crate::{Data, DataMut, DataOwned, RawData};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `modulus`.
    #[inline]
    pub fn sub<M, A>(mut self, rhs: &NttPolynomial<A, T>, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign(rhs, modulus);
        self
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self -= rhs` according to `modulus`.
    #[inline]
    pub fn sub_assign<M, A>(&mut self, rhs: &NttPolynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(a, &b)| modulus.reduce_sub_assign(a, b));
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut NttPolynomial<A, T>, modulus: M)
    where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(self.poly_length(), rhs.poly_length());
        debug_assert_eq!(self.poly_length(), result.poly_length());

        izip!(self.iter(), rhs.iter(), result.iter_mut())
            .for_each(|(&a, &b, c)| *c = modulus.reduce_sub(a, b))
    }
}
