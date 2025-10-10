use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceAdd, ReduceAddAssign};

use crate::{ArrayBase, Data, DataMut, DataOwned, PolyLength, RawData};

use super::CrtPolynomial;

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self + rhs` according to `moduli`.
    #[inline]
    pub fn add<M, A>(
        mut self,
        rhs: &CrtPolynomial<A, T>,
        moduli: &[M],
        poly_length: PolyLength,
    ) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign(rhs, moduli, poly_length);
        self
    }
}

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self += rhs` according to `moduli`.
    #[inline]
    pub fn add_assign<M, A>(
        &mut self,
        rhs: &CrtPolynomial<A, T>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(self.iter_mut(poly_length), rhs.iter(poly_length), moduli).for_each(
            |(xs, ys, modulus)| {
                ArrayBase(xs).add_assign(&ArrayBase(ys), *modulus);
            },
        );
    }
}

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self + rhs` according to `moduli`.
    #[inline]
    pub fn add_inplace<M, A>(
        &self,
        rhs: &Self,
        result: &mut CrtPolynomial<A, T>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter(poly_length),
            rhs.iter(poly_length),
            result.iter_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, modulus)| {
            ArrayBase(xs).add_inplace(&ArrayBase(ys), &mut ArrayBase(zs), *modulus);
        });
    }
}
