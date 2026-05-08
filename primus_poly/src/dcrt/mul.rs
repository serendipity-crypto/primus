use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{Data, DataMut, RawData, UnsignedInteger, izip};
use primus_reduce::ops::{ReduceMul, ReduceMulAdd, ReduceMulAssign};

use crate::ArrayBase;

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S>
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

    /// Performs `self += scalar * rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_scalar_assign<M, A>(
        &mut self,
        rhs: &DcrtPolynomial<A>,
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

    /// Performs `self * rhs` according to `moduli`.
    #[inline]
    pub fn mul<M, A>(mut self, rhs: &DcrtPolynomial<A>, poly_length: usize, moduli: &[M]) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.mul_assign(rhs, poly_length, moduli);
        self
    }

    /// Performs `self *= rhs` according to `moduli`.
    #[inline]
    pub fn mul_assign<M, A>(&mut self, rhs: &DcrtPolynomial<A>, poly_length: usize, moduli: &[M])
    where
        M: Copy + ReduceMulAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            rhs.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, &modulus)| {
            ArrayBase(xs).mul_element_wise_assign(&ArrayBase(ys), modulus)
        })
    }

    /// Inverse butterfly with a Shoup-factor polynomial.
    ///
    /// `(self, result) = (self + rhs, (self_orig - rhs) * w)`.
    ///
    /// `self` and `rhs` are expected in `[0, q)`. Both outputs are written
    /// back in `[0, q)`.
    #[inline]
    pub fn butterfly_mul_factor_inplace<A, B>(
        &mut self,
        rhs: &DcrtPolynomial<A>,
        w: &[ShoupFactor<T>],
        result: &mut DcrtPolynomial<B>,
        poly_length: usize,
        moduli: &[T],
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let lhs = self.0.as_mut_slice();
        let rhs = rhs.0.as_slice();
        let result = result.0.as_mut_slice();

        debug_assert_eq!(lhs.len(), rhs.len());
        debug_assert_eq!(lhs.len(), result.len());
        debug_assert_eq!(lhs.len(), w.len());
        debug_assert_eq!(lhs.len(), poly_length * moduli.len());

        moduli
            .iter()
            .enumerate()
            .for_each(|(modulus_index, &modulus)| {
                let offset = modulus_index * poly_length;
                (0..poly_length).for_each(|coeff_index| {
                    let index = offset + coeff_index;

                    // SAFETY: `offset < lhs.len()` and `coeff_index < poly_length`.
                    unsafe {
                        let a = lhs.get_unchecked_mut(index);
                        let s = *rhs.get_unchecked(index);
                        let w = *w.get_unchecked(index);
                        let b = result.get_unchecked_mut(index);

                        let a_orig = *a;
                        let sum_diff = modulus - s;
                        *a = if sum_diff > a_orig {
                            a_orig + s
                        } else {
                            a_orig.wrapping_sub(sum_diff)
                        };
                        let diff = if s > a_orig {
                            a_orig.wrapping_sub(s).wrapping_add(modulus)
                        } else {
                            a_orig - s
                        };
                        *b = w.factor_mul_modulo(diff, modulus);
                    }
                });
            })
    }
}

impl<S, T> DcrtPolynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs `result = self * rhs` according to `moduli`.
    #[inline]
    pub fn mul_inplace<M, A, B>(
        &self,
        rhs: &DcrtPolynomial<A>,
        result: &mut DcrtPolynomial<B>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMul<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
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

    /// Performs `result = self * scalar` according to `moduli`.
    #[inline]
    pub fn mul_factor_inplace<A>(
        &self,
        scalar: &[ShoupFactor<T>],
        result: &mut DcrtPolynomial<A>,
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
