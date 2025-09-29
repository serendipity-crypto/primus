use core::cmp::Ordering;

use crate::{UnsignedInteger, izip};

/// A trait for big integer types represented as slices of smaller unsigned integer types.
pub trait BigInteger: AsRef<[Self::ValueT]> + AsMut<[Self::ValueT]> {
    /// The underlying unsigned integer type used in the slice representation.
    type ValueT;

    /// Gets the bits count of the big integer slice.
    #[must_use]
    fn bits_count(&self) -> u32;
}

/// Implement BigInteger for slices of any UnsignedInteger type.
impl<T: UnsignedInteger> BigInteger for [T] {
    type ValueT = T;

    #[inline]
    fn bits_count(&self) -> u32 {
        self.iter()
            .enumerate()
            .rev()
            .find(|(_, v)| !v.is_zero())
            .map_or(0, |(i, v)| T::BITS * (i as u32 + 1) - v.leading_zeros())
    }
}

/// A trait providing various operations on big integers represented as slices of unsigned integers.
pub trait BigIntegerOps: BigInteger {
    /// Left shifts the big integer slice.
    fn slice_left_shift_assign(&mut self, bits: u32);

    /// Left shifts the big integer slice.
    fn slice_right_shift_assign(&mut self, bits: u32);

    /// Adds a value to the big integer slice, returning true if there was a carry.
    #[must_use]
    fn slice_add_value_assign(&mut self, value: Self::ValueT) -> bool;

    /// Subtracts a value from the big integer slice, returning true if there was a borrow.
    #[must_use]
    fn slice_sub_value_assign(&mut self, value: Self::ValueT) -> bool;

    /// Multiplies the big integer slice by a value, returning any carry that results.
    #[must_use]
    fn slice_mul_value_assign(&mut self, value: Self::ValueT) -> Self::ValueT;

    /// Adds a value to the big integer slice, returning true if there was a carry.
    #[must_use]
    fn slice_add_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool;

    /// Multiplies the big integer slice by a value, storing the result in another slice.
    #[must_use]
    fn slice_mul_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT;

    /// Multiplies the big integer slice by a value, then add to another slice.
    #[must_use]
    fn slice_mul_value_add_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT;

    /// Adds another big integer slice to this one, returning true if there was a carry.
    #[must_use]
    fn slice_add_assign(&mut self, other: &Self) -> bool;

    /// Subtracts another big integer slice from this one, returning true if there was a borrow.
    #[must_use]
    fn slice_sub_assign(&mut self, other: &Self) -> bool;

    /// Adds two big integer slices to result, returning true if there was a carry.
    #[must_use]
    fn slice_add_inplace(&self, other: &Self, result: &mut Self) -> bool;

    /// Subtracts another big integer slice from this one, returning true if there was a borrow.
    #[must_use]
    fn slice_sub_inplace(&self, other: &Self, result: &mut Self) -> bool;

    /// Compares this big integer slice with another, returning an Ordering.
    #[must_use]
    fn slice_cmp(&self, other: &Self) -> Ordering;

    /// Adds another big integer slice to this one modulo a given modulus.
    fn slice_add_modulo_assign(&mut self, other: &Self, modulus: &Self);

    /// Subs another big integer slice to this one modulo a given modulus.
    fn slice_sub_modulo_assign(&mut self, other: &Self, modulus: &Self);

    fn slice_neg_modulo_assign(&mut self, modulus: &Self);

    /// Adds two big integer slices to result modulo a given modulus.
    fn slice_add_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self);

    /// Subs another big integer slice to this one modulo a given modulus.
    fn slice_sub_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self);

    fn slice_neg_modulo_inplace(&self, result: &mut Self, modulus: &Self);
}

impl<T: UnsignedInteger> BigIntegerOps for [T] {
    #[inline]
    fn slice_left_shift_assign(&mut self, bits: u32) {
        if bits != 0 {
            let mut pre = T::ZERO;
            let mut temp = T::ZERO;
            let right_shift_bits = T::BITS - bits;
            self.iter_mut().for_each(|value| {
                temp = *value;
                *value = *value << bits | pre >> right_shift_bits;
                pre = temp;
            });
        }
    }

    #[inline]
    fn slice_right_shift_assign(&mut self, bits: u32) {
        if bits != 0 {
            let mut pre = T::ZERO;
            let mut temp = T::ZERO;
            let left_shift_bits = T::BITS - bits;
            self.iter_mut().rev().for_each(|value| {
                temp = *value;
                *value = pre << left_shift_bits | *value >> bits;
                pre = temp;
            });
        }
    }

    #[inline]
    fn slice_add_value_assign(&mut self, value: Self::ValueT) -> bool {
        let mut carry;
        match self {
            [first, other @ ..] => {
                (*first, carry) = first.overflowing_add(value);
                for v in other.iter_mut() {
                    if !carry {
                        return false;
                    }
                    (*v, carry) = (*v).overflowing_add(T::ONE);
                }
                carry
            }
            _ => unreachable!(),
        }
    }

    #[inline]
    fn slice_sub_value_assign(&mut self, value: Self::ValueT) -> bool {
        let mut borrow;
        match self {
            [first, other @ ..] => {
                (*first, borrow) = first.overflowing_sub(value);
                for v in other.iter_mut() {
                    if !borrow {
                        return false;
                    }
                    (*v, borrow) = (*v).overflowing_sub(T::ONE);
                }
                borrow
            }
            _ => unreachable!(),
        }
    }

    #[inline]
    fn slice_mul_value_assign(&mut self, value: Self::ValueT) -> Self::ValueT {
        if value.is_zero() {
            self.fill(T::ZERO);
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for ele in self.iter_mut() {
            (*ele, carry) = value.carrying_mul(*ele, carry);
        }

        carry
    }

    fn slice_add_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool {
        debug_assert_eq!(self.len(), result.len());

        let mut carry;

        let mut a_iter = self.iter();
        let mut b_iter = result.iter_mut();

        let a_first = a_iter.next().unwrap();
        let b_first = b_iter.next().unwrap();

        (*b_first, carry) = a_first.overflowing_add(value);

        while carry {
            if let Some(a_next) = a_iter.next()
                && let Some(b_next) = b_iter.next()
            {
                (*b_next, carry) = a_next.overflowing_add(T::ONE);
            } else {
                return carry;
            }
        }

        for (a, b) in a_iter.zip(b_iter) {
            *b = *a;
        }

        carry
    }

    #[inline]
    fn slice_mul_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT {
        debug_assert_eq!(result.len(), self.len());

        if value.is_zero() {
            result.fill(T::ZERO);
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for (ele, res) in self.iter().zip(result.iter_mut()) {
            (*res, carry) = value.carrying_mul(*ele, carry);
        }

        carry
    }

    #[inline]
    fn slice_mul_value_add_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT {
        debug_assert_eq!(result.len(), self.len());

        if value.is_zero() {
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for (ele, res) in self.iter().zip(result.iter_mut()) {
            (*res, carry) = value.carrying_mul_add(*ele, *res, carry);
        }

        carry
    }

    #[inline]
    fn slice_add_assign(&mut self, other: &Self) -> bool {
        debug_assert_eq!(self.len(), other.len());

        let mut carry = false;

        for (xs, ys) in self.iter_mut().zip(other) {
            (*xs, carry) = xs.carrying_add(*ys, carry);
        }

        carry
    }

    #[inline]
    fn slice_sub_assign(&mut self, other: &Self) -> bool {
        debug_assert_eq!(self.len(), other.len());

        let mut borrow = false;

        for (xs, ys) in self.iter_mut().zip(other) {
            (*xs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        borrow
    }

    #[inline]
    fn slice_add_inplace(&self, other: &Self, result: &mut Self) -> bool {
        debug_assert_eq!(self.len(), other.len());
        debug_assert_eq!(self.len(), result.len());

        let mut carry = false;
        for (xs, ys, zs) in izip!(self, other, result) {
            (*zs, carry) = xs.carrying_add(*ys, carry);
        }

        carry
    }

    #[inline]
    fn slice_sub_inplace(&self, other: &Self, result: &mut Self) -> bool {
        debug_assert_eq!(self.len(), other.len());
        debug_assert_eq!(self.len(), result.len());

        let mut borrow = false;
        for (xs, ys, zs) in izip!(self, other, result) {
            (*zs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        borrow
    }

    #[inline]
    fn slice_cmp(&self, other: &Self) -> Ordering {
        debug_assert_eq!(self.len(), other.len());

        for (a, b) in self.iter().rev().zip(other.iter().rev()) {
            match a.cmp(b) {
                Ordering::Equal => continue,
                non_eq => return non_eq,
            }
        }

        Ordering::Equal
    }

    #[inline]
    fn slice_add_modulo_assign(&mut self, other: &Self, modulus: &Self) {
        debug_assert!(self.len() == other.len() && self.len() == modulus.len());
        debug_assert!(self.slice_cmp(modulus).is_lt());
        debug_assert!(other.slice_cmp(modulus).is_lt());

        let carry = self.slice_add_assign(other);
        if carry || self.slice_cmp(modulus).is_ge() {
            let _ = self.slice_sub_assign(modulus);
        }
    }

    #[inline]
    fn slice_sub_modulo_assign(&mut self, other: &Self, modulus: &Self) {
        debug_assert!(self.len() == other.len() && self.len() == modulus.len());
        debug_assert!(self.slice_cmp(modulus).is_lt());
        debug_assert!(other.slice_cmp(modulus).is_lt());

        if self.slice_sub_assign(other) {
            let _ = self.slice_add_assign(modulus);
        }
    }

    #[inline]
    fn slice_neg_modulo_assign(&mut self, modulus: &Self) {
        debug_assert!(self.len() == modulus.len());
        debug_assert!(self.slice_cmp(modulus).is_lt());

        if !self.iter().all(T::is_zero) {
            let mut borrow = false;
            for (xs, ys) in self.iter_mut().zip(modulus) {
                (*xs, borrow) = ys.borrowing_sub(*xs, borrow);
            }
        }
    }

    #[inline]
    fn slice_add_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self) {
        debug_assert!(
            self.len() == other.len() && self.len() == result.len() && self.len() == modulus.len()
        );
        debug_assert!(self.slice_cmp(modulus).is_lt());
        debug_assert!(other.slice_cmp(modulus).is_lt());

        let carry = self.slice_add_inplace(other, result);
        if carry || result.slice_cmp(modulus).is_ge() {
            let _ = result.slice_sub_assign(modulus);
        }
    }

    #[inline]
    fn slice_sub_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self) {
        debug_assert!(
            self.len() == other.len() && self.len() == result.len() && self.len() == modulus.len()
        );
        debug_assert!(self.slice_cmp(modulus).is_lt());
        debug_assert!(other.slice_cmp(modulus).is_lt());

        if self.slice_sub_inplace(other, result) {
            let _ = result.slice_add_assign(modulus);
        }
    }

    #[inline]
    fn slice_neg_modulo_inplace(&self, result: &mut Self, modulus: &Self) {
        debug_assert!(self.len() == result.len() && self.len() == modulus.len());
        debug_assert!(self.slice_cmp(modulus).is_lt());

        if self.iter().all(T::is_zero) {
            result.fill(T::ZERO);
        } else {
            let mut borrow = false;
            for (xs, ys, zs) in izip!(self, modulus, result) {
                (*zs, borrow) = ys.borrowing_sub(*xs, borrow);
            }
        }
    }
}

/// Multiplies many values together, returning the result as a big integer slice.
pub fn multiply_many_values<T: UnsignedInteger>(values: &[T]) -> Vec<T> {
    let mut result = Vec::with_capacity(values.len());
    result.push(values[0]);
    for &v in values.iter().skip(1) {
        let carry = result.slice_mul_value_assign(v);
        if !carry.is_zero() {
            result.push(carry);
        }
    }
    result.shrink_to_fit();
    result
}

/// Multiplies many values together, except for one specified by index, returning the result as a big integer slice.
pub fn multiply_many_values_except<T: UnsignedInteger>(values: &[T], except: usize) -> Vec<T> {
    let mut result = Vec::with_capacity(values.len() - 1);
    result.push(T::ONE);

    for (i, &v) in values.iter().enumerate() {
        if i == except {
            continue;
        }
        let carry = result.as_mut_slice().slice_mul_value_assign(v);
        if !carry.is_zero() {
            result.push(carry);
        }
    }

    result.shrink_to_fit();
    result
}

/// Multiplies many values together, except for one specified by index, returning the result as a big integer slice.
pub fn multiply_many_values_except_inplace<T: UnsignedInteger>(
    values: &[T],
    except: usize,
    result: &mut [T],
) {
    result.fill(T::ZERO);
    result[0] = T::ONE;
    let mut len = 1;

    for (_, &v) in values.iter().enumerate().filter(|(i, _)| *i != except) {
        let carry = result[0..len].slice_mul_value_assign(v);
        if !carry.is_zero() {
            result[len] = carry;
            len += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use rand::{
        Rng,
        distr::{Distribution, Uniform},
    };

    use super::*;

    type ValueT = u32;

    fn compose(value: &[ValueT]) -> u128 {
        assert!(value.len() <= 4);
        let mut result = 0u128;
        for &r in value.into_iter().rev() {
            result <<= ValueT::BITS;
            result |= r as u128;
        }
        result
    }

    #[test]
    fn test_big_integer_ops() {
        let mut rng = rand::rng();
        let moduli: [ValueT; 3] = [134215681, 134176769, 132120577];
        let composed_modulus = multiply_many_values(&moduli);
        let m = compose(&composed_modulus);

        let bits_count = composed_modulus.bits_count();
        assert_eq!(bits_count, 128 - m.leading_zeros());

        let distrs = moduli.map(|m| Uniform::new(0, m).unwrap());

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a_vec = multiply_many_values(&a_residues);
        let mut a = compose(&a_vec);

        a_vec.slice_right_shift_assign(3);
        a >>= 3;
        assert_eq!(a, compose(&a_vec));

        a_vec.slice_left_shift_assign(3);
        a <<= 3;
        assert_eq!(a, compose(&a_vec));

        let v: ValueT = rng.random();
        let _r = a_vec.slice_add_value_assign(v);
        a += v as u128;
        assert_eq!(a, compose(&a_vec));

        let _r = a_vec.slice_sub_value_assign(v);
        a -= v as u128;
        assert_eq!(a, compose(&a_vec));

        let r = a_vec.slice_mul_value_assign(v);
        let mut a_vecp = a_vec.clone();
        a_vecp.push(r);
        a *= v as u128;
        assert_eq!(a, compose(&a_vecp));

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let b_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a_vec = multiply_many_values(&a_residues);
        let b_vec = multiply_many_values(&b_residues);
        let a = compose(&a_vec);
        let b = compose(&b_vec);

        let _r = a_vec.slice_add_assign(&b_vec);
        assert_eq!(a + b, compose(&a_vec));

        let _r = a_vec.slice_sub_assign(&b_vec);
        assert_eq!(a, compose(&a_vec));

        a_vec.slice_add_modulo_assign(&b_vec, &composed_modulus);
        let r = (a + b) % m;
        assert_eq!(r, compose(&a_vec));

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let b_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a_vec = multiply_many_values(&a_residues);
        let b_vec = multiply_many_values(&b_residues);
        let a = compose(&a_vec);
        let b = compose(&b_vec);

        a_vec.slice_sub_modulo_assign(&b_vec, &composed_modulus);
        let r = (a + m - b) % m;
        assert_eq!(r, compose(&a_vec))
    }
}
