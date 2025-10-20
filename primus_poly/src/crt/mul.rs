use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceMulAdd, ReduceMulAssign};

use crate::{ArrayBase, Data, DataMut, RawData};

use super::CrtPolynomial;

impl<S, T> CrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, poly_length: usize, moduli: &[M]) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, poly_length, moduli);
        self
    }

    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_each_modulus_mut(poly_length)
            .zip(moduli)
            .for_each(|(poly, &modulus)| ArrayBase(poly).mul_scalar_assign(scalar, modulus))
    }

    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_residues_assign<M>(
        &mut self,
        scalar_residues: &[T],
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMulAssign<T>,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            scalar_residues,
            moduli
        )
        .for_each(|(poly, &scalar, &modulus)| ArrayBase(poly).mul_scalar_assign(scalar, modulus))
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M, A>(
        &mut self,
        rhs: &CrtPolynomial<A, T>,
        scalar: T,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, &modulus)| {
            ArrayBase(xs).add_mul_scalar_assign(&ArrayBase(ys), scalar, modulus);
        });
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_residues_assign<M, A>(
        &mut self,
        rhs: &CrtPolynomial<A, T>,
        scalar_residues: &[T],
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            scalar_residues,
            moduli
        )
        .for_each(|(xs, ys, &scalar, &modulus)| {
            ArrayBase(xs).add_mul_scalar_assign(&ArrayBase(ys), scalar, modulus);
        });
    }
}
