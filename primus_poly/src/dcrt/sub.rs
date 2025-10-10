use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceSub, ReduceSubAssign};

use crate::{ArrayBase, Data, DataMut, DataOwned, PolyLength, RawData};

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self - rhs` according to `moduli`.
    #[inline]
    pub fn sub<M, A>(mut self, rhs: &Self, moduli: &[M], poly_length: PolyLength) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign(rhs, moduli, poly_length);
        self
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self -= rhs` according to `moduli`.
    #[inline]
    pub fn sub_assign<M, A>(
        &mut self,
        rhs: &DcrtPolynomial<A, T>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(self.iter_mut(poly_length), rhs.iter(poly_length), moduli).for_each(
            |(xs, ys, modulus)| {
                ArrayBase(xs).sub_assign(&ArrayBase(ys), *modulus);
            },
        );
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self - rhs` according to `moduli`.
    #[inline]
    pub fn sub_inplace<M, A>(
        &self,
        rhs: &Self,
        result: &mut DcrtPolynomial<A, T>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter(poly_length),
            rhs.iter(poly_length),
            result.iter_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, modulus)| {
            ArrayBase(xs).sub_inplace(&ArrayBase(ys), &mut ArrayBase(zs), *modulus);
        });
    }
}
