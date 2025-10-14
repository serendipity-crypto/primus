use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceMul, ReduceMulAdd, ReduceMulAssign};

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Performs `self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, moduli: &[M], poly_length: usize) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, moduli, poly_length);
        self
    }

    /// Performs `self * rhs` according to `moduli`.
    #[inline]
    pub fn mul<M, A>(mut self, rhs: &DcrtPolynomial<A, T>, moduli: &[M], poly_length: usize) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.mul_assign(rhs, moduli, poly_length);
        self
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, moduli: &[M], poly_length: usize)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_each_modulus_mut(poly_length)
            .zip(moduli)
            .for_each(|(poly, modulus)| ArrayBase(poly).mul_scalar_assign(scalar, *modulus))
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M>(
        &mut self,
        rhs: &Self,
        scalar: T,
        moduli: &[M],
        poly_length: usize,
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, modulus)| {
            ArrayBase(xs).add_mul_scalar_assign(&ArrayBase(ys), scalar, *modulus);
        });
    }

    /// Performs `self *= rhs` according to `moduli`.
    #[inline]
    pub fn mul_assign<M, A>(&mut self, rhs: &DcrtPolynomial<A, T>, moduli: &[M], poly_length: usize)
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, modulus)| {
            ArrayBase(xs).mul_element_wise_assign(&ArrayBase(ys), *modulus)
        })
    }

    /// Performs `result = self * rhs` according to `moduli`.
    #[inline]
    pub fn mul_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M], poly_length: usize)
    where
        M: Copy + ReduceMul<T, Output = T>,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            rhs.iter_each_modulus(poly_length),
            result.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, modulus)| {
            ArrayBase(xs).mul_element_wise_inplace(&ArrayBase(ys), &mut ArrayBase(zs), *modulus);
        })
    }
}
