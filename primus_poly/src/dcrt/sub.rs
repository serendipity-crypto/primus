use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceSub, ReduceSubAssign};

use crate::{ArrayBase, Data, DataMut, RawData};

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `moduli`.
    #[inline]
    pub fn sub<M, A>(mut self, rhs: &Self, poly_length: usize, moduli: &[M]) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign(rhs, poly_length, moduli);
        self
    }

    /// Performs `self -= rhs` according to `moduli`.
    #[inline]
    pub fn sub_assign<M, A>(&mut self, rhs: &DcrtPolynomial<A, T>, poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, &modulus)| {
            ArrayBase(xs).sub_element_wise_assign(&ArrayBase(ys), modulus);
        });
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_inplace<M, A, B>(
        &self,
        rhs: &DcrtPolynomial<A, T>,
        result: &mut DcrtPolynomial<B, T>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            rhs.iter_each_modulus(poly_length),
            result.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, &modulus)| {
            ArrayBase(xs).sub_element_wise_inplace(&ArrayBase(ys), &mut ArrayBase(zs), modulus);
        });
    }

    /// Performs `rhs = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_to_right<M, A>(
        &self,
        rhs: &mut DcrtPolynomial<A, T>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            rhs.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, &modulus)| {
            ArrayBase(xs).sub_element_wise_to_right(&mut ArrayBase(ys), modulus);
        });
    }
}
