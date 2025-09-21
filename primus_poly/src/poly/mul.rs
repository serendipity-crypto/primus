use integer::UnsignedInteger;
use num_traits::{ConstZero, One, Zero};
use primus_factor::{FactorMul, ShoupFactor};
use reduce::ops::*;
use uint_modulus::UintModulus;

use super::Polynomial;

impl<T> Polynomial<T>
where
    T: Copy + Zero + ConstZero + One + PartialEq,
{
    /// Multiply `self` with a scalar.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Multiply `self` with a scalar assign.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        if scalar.is_zero() {
            self.set_zero();
        } else if !scalar.is_one() {
            self.iter_mut()
                .for_each(|v| modulus.reduce_mul_assign(v, scalar))
        }
    }

    /// Add the multiply result `rhs` with a scalar inplace.
    #[inline]
    pub fn add_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, modulus: M)
    where
        M: Copy + ReduceAddAssign<T> + ReduceMulAdd<T, Output = T>,
    {
        if scalar.is_one() {
            self.add_assign(rhs, modulus);
        } else if !scalar.is_zero() {
            self.iter_mut()
                .zip(rhs.iter())
                .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
        }
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

        let coeff_count = self.coeff_count();
        debug_assert_eq!(coeff_count, poly2.len());
        debug_assert_eq!(coeff_count, destination.coeff_count());

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
    /// Multiply `self` with a scalar.
    #[inline]
    pub fn mul_factor_scalar(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_scalar_assign(scalar, modulus);
        self
    }

    /// Multiply `self` with a scalar inplace.
    #[inline]
    pub fn mul_factor_scalar_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        if scalar.value().is_zero() {
            self.set_zero();
        } else if !scalar.value().is_one() {
            self.iter_mut()
                .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
        }
    }

    /// Multiply `self` with the a shoup scalar and add to self.
    #[inline]
    pub fn add_mul_factor_scalar_assign(&mut self, rhs: &Self, scalar: ShoupFactor<T>, modulus: T) {
        if scalar.value().is_one() {
            self.add_assign(rhs, UintModulus(modulus));
        } else if !scalar.value().is_zero() {
            self.iter_mut().zip(rhs).for_each(|(r, &v)| {
                UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
            });
        }
    }
}
