use primus_integer::UnsignedInteger;
use primus_reduce::{lazy_ops::*, ops::*};

use crate::UintModulus;

use super::MontgomeryModulus;

// LazyReduce implementations
impl<T: UnsignedInteger> LazyReduce<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn lazy_reduce(self, value: T) -> T {
        let (lo, hi) = value.widening_mul(self.r);
        self.montgomery_reduce([lo, hi])
    }
}

impl<T: UnsignedInteger> LazyReduceAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn lazy_reduce_assign(self, value: &mut T) {
        *value = self.lazy_reduce(*value);
    }
}

impl<T: UnsignedInteger> LazyReduceMul<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn lazy_reduce_mul(self, a: T, b: T) -> Self::Output {
        let (lo, hi) = a.widening_mul(b);
        self.montgomery_reduce([lo, hi])
    }
}

impl<T: UnsignedInteger> LazyReduceMulAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn lazy_reduce_mul_assign(self, a: &mut T, b: T) {
        let (lo, hi) = a.widening_mul(b);
        *a = self.montgomery_reduce([lo, hi]);
    }
}

impl<T: UnsignedInteger> LazyReduceMulAdd<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn lazy_reduce_mul_add(self, a: T, b: T, c: T) -> Self::Output {
        self.reduce_add(self.reduce_mul(a, b), c)
    }
}

impl<T: UnsignedInteger> LazyReduceMulAddAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn lazy_reduce_mul_add_assign(self, a: &mut T, b: T, c: T) {
        *a = self.lazy_reduce_mul_add(*a, b, c);
    }
}

// Reduce implementations
impl<T: UnsignedInteger> Reduce<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce(self, value: T) -> Self::Output {
        let (lo, hi) = value.widening_mul(self.r);
        self.montgomery_reduce([lo, hi])
    }
}

impl<T: UnsignedInteger> ReduceAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_assign(self, value: &mut T) {
        *value = self.reduce(*value);
    }
}

impl<T: UnsignedInteger> ReduceOnce<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_once(self, value: T) -> Self::Output {
        UintModulus(self.value).reduce_once(value)
    }
}

impl<T: UnsignedInteger> ReduceOnceAssign<T> for MontgomeryModulus<T> {
    #[inline(always)]
    fn reduce_once_assign(self, value: &mut T) {
        UintModulus(self.value).reduce_once_assign(value);
    }
}

impl<T: UnsignedInteger> ReduceAdd<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_add(self, a: T, b: T) -> Self::Output {
        UintModulus(self.value).reduce_add(a, b)
    }
}

impl<T: UnsignedInteger> ReduceAddAssign<T> for MontgomeryModulus<T> {
    #[inline(always)]
    fn reduce_add_assign(self, a: &mut T, b: T) {
        UintModulus(self.value).reduce_add_assign(a, b);
    }
}

impl<T: UnsignedInteger> ReduceDouble<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_double(self, value: T) -> Self::Output {
        UintModulus(self.value).reduce_double(value)
    }
}

impl<T: UnsignedInteger> ReduceDoubleAssign<T> for MontgomeryModulus<T> {
    #[inline(always)]
    fn reduce_double_assign(self, value: &mut T) {
        UintModulus(self.value).reduce_double_assign(value);
    }
}

impl<T: UnsignedInteger> ReduceSub<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_sub(self, a: T, b: T) -> Self::Output {
        UintModulus(self.value).reduce_sub(a, b)
    }
}

impl<T: UnsignedInteger> ReduceSubAssign<T> for MontgomeryModulus<T> {
    #[inline(always)]
    fn reduce_sub_assign(self, a: &mut T, b: T) {
        UintModulus(self.value).reduce_sub_assign(a, b);
    }
}

impl<T: UnsignedInteger> ReduceNeg<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline(always)]
    fn reduce_neg(self, value: T) -> Self::Output {
        UintModulus(self.value).reduce_neg(value)
    }
}

impl<T: UnsignedInteger> ReduceNegAssign<T> for MontgomeryModulus<T> {
    #[inline(always)]
    fn reduce_neg_assign(self, value: &mut T) {
        UintModulus(self.value).reduce_neg_assign(value);
    }
}

impl<T: UnsignedInteger> ReduceMul<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn reduce_mul(self, a: T, b: T) -> Self::Output {
        let (lo, hi) = a.widening_mul(b);
        self.montgomery_reduce([lo, hi])
    }
}

impl<T: UnsignedInteger> ReduceMulAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_mul_assign(self, a: &mut T, b: T) {
        let (lo, hi) = a.widening_mul(b);
        *a = self.montgomery_reduce([lo, hi]);
    }
}

impl<T: UnsignedInteger> ReduceSquare<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn reduce_square(self, value: T) -> Self::Output {
        let (lo, hi) = value.widening_mul(value);
        self.montgomery_reduce([lo, hi])
    }
}

impl<T: UnsignedInteger> ReduceSquareAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_square_assign(self, value: &mut T) {
        let (lo, hi) = value.widening_mul(*value);
        *value = self.montgomery_reduce([lo, hi]);
    }
}

impl<T: UnsignedInteger> ReduceMulAdd<T> for MontgomeryModulus<T> {
    type Output = T;

    #[inline]
    fn reduce_mul_add(self, a: T, b: T, c: T) -> Self::Output {
        self.reduce_add(self.reduce_mul(a, b), c)
    }
}

impl<T: UnsignedInteger> ReduceMulAddAssign<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_mul_add_assign(self, a: &mut T, b: T, c: T) {
        *a = self.reduce_mul_add(*a, b, c);
    }
}

impl<T: UnsignedInteger> ReduceExp<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_exp<E: UnsignedInteger>(self, base: T, mut exp: E) -> T {
        if exp.is_zero() {
            return T::ONE;
        }

        if base.is_zero() {
            return T::ZERO;
        }

        debug_assert!(base < self.value);

        let mut power: T = base;

        let exp_trailing_zeros = exp.trailing_zeros();
        if exp_trailing_zeros > 0 {
            for _ in 0..exp_trailing_zeros {
                self.reduce_square_assign(&mut power);
            }
            exp >>= exp_trailing_zeros;
        }

        if exp.is_one() {
            return power;
        }

        let mut intermediate: T = power;
        for _ in 1..(E::BITS - exp.leading_zeros()) {
            exp >>= 1;
            self.reduce_square_assign(&mut power);
            if !(exp & E::ONE).is_zero() {
                self.reduce_mul_assign(&mut intermediate, power);
            }
        }
        intermediate
    }
}

impl<T: UnsignedInteger> ReduceExpPowOf2<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_exp_power_of_2(self, base: T, exp_log: u32) -> T {
        if base.is_zero() {
            return T::ZERO;
        }

        let mut power = base;

        for _ in 0..exp_log {
            self.reduce_square_assign(&mut power);
        }

        power
    }
}

impl<T: UnsignedInteger> ReduceDotProduct<T> for MontgomeryModulus<T> {
    #[inline]
    fn reduce_dot_product(self, a: impl AsRef<[T]>, b: impl AsRef<[T]>) -> T {
        a.as_ref()
            .iter()
            .zip(b.as_ref())
            .fold(T::ZERO, |acc, (&x, &y)| self.reduce_mul_add(x, y, acc))
    }

    #[inline]
    fn reduce_dot_product_iter(
        self,
        a: impl IntoIterator<Item = T>,
        b: impl IntoIterator<Item = T>,
    ) -> T {
        a.into_iter()
            .zip(b)
            .fold(T::ZERO, |acc, (x, y)| self.reduce_mul_add(x, y, acc))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_mul() {
        let m = MontgomeryModulus::<u32>::new(17);
        let a = 10u32;
        let b = 15u32;

        let result = m.reduce_mul(a, b);
        assert_eq!(result, (10 * 15) % 17);
    }

    #[test]
    fn test_reduce_square() {
        let m = MontgomeryModulus::<u32>::new(17);
        let value = 10u32;

        let result = m.reduce_square(value);
        assert_eq!(result, (10 * 10) % 17);
    }

    #[test]
    fn test_reduce_exp() {
        let m = MontgomeryModulus::<u32>::new(17);
        let base = 3u32;
        let exp = 5u32;

        let result = m.reduce_exp(base, exp);

        // 3^5 = 243, 243 % 17 = 5
        assert_eq!(result, 5);
    }

    #[test]
    fn test_montgomery_reduction() {
        let m = MontgomeryModulus::<u32>::new(17);

        // Test that montgomery_reduce works correctly
        let value = [100u32, 0u32];
        let result = m.montgomery_reduce(value);
        assert!(result < 17);
    }
}
