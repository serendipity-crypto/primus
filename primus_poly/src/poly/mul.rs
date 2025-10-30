use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::UnsignedInteger;
use primus_modulus::UintModulus;
use primus_reduce::ops::*;

use crate::{Data, DataMut, RawData};

use super::Polynomial;

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_scalar<M>(mut self, scalar: T, modulus: M) -> Self
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.mul_scalar_assign(scalar, modulus);
        self
    }

    /// Performs `self * scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor(mut self, scalar: ShoupFactor<T>, modulus: T) -> Self {
        self.mul_factor_assign(scalar, modulus);
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
    pub fn add_mul_scalar_assign<M, A>(&mut self, rhs: &Polynomial<A, T>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut()
            .zip(rhs.iter())
            .for_each(|(r, &v)| *r = modulus.reduce_mul_add(v, scalar, *r));
    }

    /// Performs `self *= scalar` according to `modulus`.
    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        self.iter_mut()
            .for_each(|v| *v = scalar.factor_mul_modulo(*v, modulus))
    }

    /// Performs `self += scalar * rhs` according to `modulus`.
    #[inline]
    pub fn add_mul_factor_assign<A>(
        &mut self,
        rhs: &Polynomial<A, T>,
        scalar: ShoupFactor<T>,
        modulus: T,
    ) where
        A: RawData<Elem = T> + Data,
    {
        self.iter_mut().zip(rhs.iter()).for_each(|(r, &v)| {
            UintModulus(modulus).reduce_add_assign(r, scalar.factor_mul_modulo(v, modulus))
        });
    }
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// A naive multiplication over polynomial.
    pub fn naive_mul_inplace<M, A, B>(
        &self,
        rhs: &Polynomial<A, T>,
        result: &mut Polynomial<B, T>,
        modulus: M,
    ) where
        M: Copy + ReduceSubAssign<T> + ReduceMul<T, Output = T> + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let a: &[T] = self.as_ref();
        let b: &[T] = rhs.as_ref();
        let c: &mut [T] = result.as_mut();

        let coeff_count = a.len();
        debug_assert_eq!(coeff_count, b.len());
        debug_assert_eq!(coeff_count, c.len());

        for i in 0..coeff_count {
            for j in 0..=i {
                c[i] = modulus.reduce_mul_add(a[j], b[i - j], c[i]);
            }
        }

        // mod (x^n + 1)
        for i in coeff_count..coeff_count * 2 - 1 {
            let k = i - coeff_count;
            for j in i - coeff_count + 1..coeff_count {
                modulus.reduce_sub_assign(&mut c[k], modulus.reduce_mul(a[j], b[i - j]));
            }
        }
    }
}
