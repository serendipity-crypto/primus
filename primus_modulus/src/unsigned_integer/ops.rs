use primus_factor::{FactorMul, LazyFactorMul};
use primus_gcd::Xgcd;
use primus_integer::UnsignedInteger;
use primus_reduce::{ReduceError, lazy_ops::LazyReduceMul, ops::*};

use super::UintModulus;

impl<T: UnsignedInteger> ReduceOnce<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_once(self, value: T) -> Self::Output {
        if value >= self.0 {
            value - self.0
        } else {
            value
        }
    }
}

impl<T: UnsignedInteger> ReduceOnceAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_once_assign(self, value: &mut T) {
        if *value >= self.0 {
            *value -= self.0;
        }
    }
}

impl<T: UnsignedInteger> ReduceAdd<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_add(self, a: T, b: T) -> Self::Output {
        let diff = self.0 - b;
        if diff > a {
            a + b
        } else {
            a.wrapping_sub(diff)
        }
    }
}

impl<T: UnsignedInteger> ReduceAddAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_add_assign(self, a: &mut T, b: T) {
        let diff = self.0 - b;
        if diff > *a {
            *a += b;
        } else {
            *a = a.wrapping_sub(diff)
        }
    }
}

impl<T: UnsignedInteger> ReduceDouble<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_double(self, value: T) -> Self::Output {
        self.reduce_add(value, value)
    }
}

impl<T: UnsignedInteger> ReduceDoubleAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_double_assign(self, value: &mut T) {
        self.reduce_add_assign(value, *value);
    }
}

impl<T: UnsignedInteger> ReduceSub<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_sub(self, a: T, b: T) -> Self::Output {
        if b > a {
            a.wrapping_sub(b).wrapping_add(self.0)
        } else {
            a - b
        }
    }
}

impl<T: UnsignedInteger> ReduceSubAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_sub_assign(self, a: &mut T, b: T) {
        if b > *a {
            *a = a.wrapping_sub(b).wrapping_add(self.0);
        } else {
            *a -= b;
        }
    }
}

impl<T: UnsignedInteger> ReduceNeg<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_neg(self, value: T) -> Self::Output {
        if value.is_zero() {
            T::ZERO
        } else {
            self.0 - value
        }
    }
}

impl<T: UnsignedInteger> ReduceNegAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_neg_assign(self, value: &mut T) {
        if !value.is_zero() {
            *value = self.0 - *value;
        }
    }
}

impl<T: UnsignedInteger> ReduceInv<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_inv(self, value: T) -> Self::Output {
        debug_assert!(self.0 > value);

        let (inv, gcd) = Xgcd::gcdinv(value, self.0);
        assert_eq!(gcd, T::ONE, "No {value}^(-1) mod {}", self.0);

        inv
    }
}

impl<T: UnsignedInteger> ReduceInvAssign<T> for UintModulus<T> {
    #[inline(always)]
    fn reduce_inv_assign(self, value: &mut T) {
        *value = self.reduce_inv(*value);
    }
}

impl<T: UnsignedInteger> TryReduceInv<T> for UintModulus<T> {
    type Output = T;

    #[inline(always)]
    fn try_reduce_inv(self, value: T) -> Result<Self::Output, ReduceError<T>> {
        debug_assert!(self.0 > value);

        let (inv, gcd) = Xgcd::gcdinv(value, self.0);

        if gcd.is_one() {
            Ok(inv)
        } else {
            Err(ReduceError::NoInverse {
                value,
                modulus: self.0,
            })
        }
    }
}

impl<T: UnsignedInteger, F> LazyReduceMul<T, F> for UintModulus<T>
where
    F: LazyFactorMul<T>,
{
    type Output = T;

    #[inline(always)]
    fn lazy_reduce_mul(self, a: T, b: F) -> Self::Output {
        b.lazy_factor_mul_modulo(a, self.0)
    }
}

impl<T: UnsignedInteger, F> ReduceMul<T, F> for UintModulus<T>
where
    F: FactorMul<T>,
{
    type Output = T;

    #[inline(always)]
    fn reduce_mul(self, a: T, b: F) -> Self::Output {
        b.factor_mul_modulo(a, self.0)
    }
}
