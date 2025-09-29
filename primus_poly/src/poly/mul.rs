use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::UnsignedInteger;
use primus_reduce::ops::*;
use primus_uint_modulus::UintModulus;

use super::Polynomial;

impl<T: Copy> Polynomial<T> {
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.iter_mut()
            .for_each(|v| modulus.reduce_mul_assign(v, scalar))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.iter_mut()
            .zip(rhs)
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }
}

impl<T: UnsignedInteger> Polynomial<T> {
    /// A naive multiplication over polynomial.
    pub fn naive_mul_inplace<M>(&self, rhs: impl AsRef<[T]>, modulus: M, destination: &mut Self)
    where
        M: Copy + ReduceSubAssign<T> + ReduceMul<T, Output = T> + ReduceMulAdd<T, Output = T>,
    {
        let poly1: &[T] = self.as_ref();
        let poly2: &[T] = rhs.as_ref();

        let coeff_count = self.poly_length();
        debug_assert_eq!(coeff_count, poly2.len());
        debug_assert_eq!(coeff_count, destination.poly_length());

        for i in 0..coeff_count {
            for j in 0..=i {
                destination[i] = modulus.reduce_mul_add(poly1[j], poly2[i - j], destination[i]);
            }
        }

        // mod (x^n + 1)
        for i in coeff_count..coeff_count * 2 - 1 {
            let k = i - coeff_count;
            for j in i - coeff_count + 1..coeff_count {
                modulus.reduce_sub_assign(
                    &mut destination[k],
                    modulus.reduce_mul(poly1[j], poly2[i - j]),
                );
            }
        }
    }
}

impl<T: UnsignedInteger> Polynomial<T> {
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_assign(scalar, modulus);
        self
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut()
            .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_factor_assign(&mut self, rhs: &Self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut().zip(rhs).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        });
    }
}
