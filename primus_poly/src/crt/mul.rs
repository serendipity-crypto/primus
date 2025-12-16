use primus_factor::ShoupFactor;
use primus_integer::{Data, DataMut, RawData, UnsignedInteger, izip};
use primus_reduce::ops::*;

use crate::ArrayBase;

use super::CrtPolynomial;

impl<S, T> CrtPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: &[T], poly_length: usize, moduli: &[M]) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, poly_length, moduli);
        self
    }

    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: &[T], poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceMulAssign<T>,
    {
        izip!(self.iter_each_modulus_mut(poly_length), scalar, moduli).for_each(
            |(poly, &scalar, &modulus)| ArrayBase(poly).mul_scalar_assign(scalar, modulus),
        )
    }

    /// Performs `self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_factor(
        mut self,
        scalar: &[ShoupFactor<T>],
        poly_length: usize,
        moduli: &[T],
    ) -> Self {
        self.mul_factor_assign(scalar, poly_length, moduli);
        self
    }

    /// Performs `self *= scalar` according to `moduli`.
    #[inline]
    pub fn mul_factor_assign(
        &mut self,
        scalar: &[ShoupFactor<T>],
        poly_length: usize,
        moduli: &[T],
    ) {
        izip!(self.iter_each_modulus_mut(poly_length), scalar, moduli).for_each(
            |(poly, &scalar, &modulus)| ArrayBase(poly).mul_factor_assign(scalar, modulus),
        )
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M, A>(
        &mut self,
        rhs: &CrtPolynomial<A>,
        scalar: &[T],
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            scalar,
            moduli
        )
        .for_each(|(xs, ys, &scalar, &modulus)| {
            ArrayBase(xs).add_mul_scalar_assign(&ArrayBase(ys), scalar, modulus);
        });
    }

    pub fn mul_monomial_assign<M>(&mut self, r: usize, poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceNegAssign<T>,
    {
        if r < poly_length {
            let rotate = |poly: &mut [T], modulus: M| {
                poly.rotate_right(r);
                poly[0..r]
                    .iter_mut()
                    .for_each(|v| modulus.reduce_neg_assign(v));
            };

            self.iter_each_modulus_mut(poly_length)
                .zip(moduli)
                .for_each(|(poly, &modulus)| rotate(poly, modulus));
        } else {
            debug_assert!(r < poly_length * 2);
            let r = r - poly_length;

            let rotate = |poly: &mut [T], modulus: M| {
                poly.rotate_right(r);
                poly[r..]
                    .iter_mut()
                    .for_each(|v| modulus.reduce_neg_assign(v));
            };

            self.iter_each_modulus_mut(poly_length)
                .zip(moduli)
                .for_each(|(poly, &modulus)| rotate(poly, modulus));
        }
    }
}

impl<S, T> CrtPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_scalar_inplace<M, A>(
        &self,
        scalar: &[T],
        result: &mut CrtPolynomial<A>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            scalar,
            result.iter_each_modulus_mut(poly_length),
            moduli
        )
        .for_each(|(in_poly, &scalar, out_poly, &modulus)| {
            ArrayBase(in_poly).mul_scalar_inplace(scalar, &mut ArrayBase(out_poly), modulus)
        })
    }

    /// Performs `result = self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_factor_inplace<A>(
        &self,
        scalar: &[ShoupFactor<T>],
        result: &mut CrtPolynomial<A>,
        poly_length: usize,
        moduli: &[T],
    ) where
        A: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.iter_each_modulus(poly_length),
            result.iter_each_modulus_mut(poly_length),
            scalar,
            moduli
        )
        .for_each(|(in_poly, out_poly, &scalar, &modulus)| {
            ArrayBase(in_poly).mul_factor_inplace(scalar, &mut ArrayBase(out_poly), modulus)
        })
    }
}
