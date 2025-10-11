use primus_integer::{UnsignedInteger, izip};
use primus_reduce::ops::{ReduceMulAdd, ReduceMulAssign};

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

use super::CrtPolynomial;

impl<S, T> CrtPolynomial<S, T>
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
}

impl<S, T> CrtPolynomial<S, T>
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
        self.0
            .chunks_exact_mut(poly_length)
            .zip(moduli)
            .for_each(|(poly, modulus)| ArrayBase(poly).mul_scalar_assign(scalar, *modulus))
    }

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M, A>(
        &mut self,
        rhs: &CrtPolynomial<A, T>,
        scalar: T,
        moduli: &[M],
        poly_length: usize,
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.0.chunks_exact_mut(poly_length),
            rhs.0.chunks_exact(poly_length),
            moduli
        )
        .for_each(|(xs, ys, modulus)| {
            ArrayBase(xs).add_mul_scalar_assign(&ArrayBase(ys), scalar, *modulus);
        });
    }
}
