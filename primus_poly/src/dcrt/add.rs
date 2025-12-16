use primus_integer::{Data, DataMut, RawData, UnsignedInteger, izip};
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use crate::ArrayBase;

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `moduli`.
    #[inline]
    pub fn add<M, A>(mut self, rhs: &DcrtPolynomial<A>, poly_length: usize, moduli: &[M]) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign(rhs, poly_length, moduli);
        self
    }

    /// Performs `self += rhs` according to `moduli`.
    #[inline]
    pub fn add_assign<M, A>(&mut self, rhs: &DcrtPolynomial<A>, poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, &modulus)| {
            ArrayBase(xs).add_element_wise_assign(&ArrayBase(ys), modulus);
        });
    }
}

impl<S, T> DcrtPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `moduli`.
    #[inline]
    pub fn add_inplace<M, A>(
        &self,
        rhs: &Self,
        result: &mut DcrtPolynomial<A>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            rhs.iter_each_modulus(poly_length),
            result.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, &modulus)| {
            ArrayBase(xs).add_element_wise_inplace(&ArrayBase(ys), &mut ArrayBase(zs), modulus);
        });
    }
}
