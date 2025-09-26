use core::cmp::Ordering;

use crate::UnsignedInteger;

/// A trait for big integer types represented as slices of smaller unsigned integer types.
pub trait BigInteger: AsRef<[Self::ValueT]> + AsMut<[Self::ValueT]> {
    /// The underlying unsigned integer type used in the slice representation.
    type ValueT;
}

/// Implement BigInteger for slices of any UnsignedInteger type.
impl<T: UnsignedInteger> BigInteger for [T] {
    type ValueT = T;
}

/// A trait providing various operations on big integers represented as slices of unsigned integers.
pub trait BigIntegerOps: BigInteger {
    /// Gets the bits count of the big integer slice.
    #[must_use]
    fn slice_value_bits_count(&self) -> u32;

    /// Left shifts the big integer slice.
    fn slice_left_shift_assign(&mut self, bits: u32);

    /// Left shifts the big integer slice.
    fn slice_right_shift_assign(&mut self, bits: u32);

    /// Adds a value to the big integer slice, returning true if there was a carry.
    #[must_use]
    fn slice_add_value_assign(&mut self, value: Self::ValueT) -> bool;

    /// Adds a value to the big integer slice, returning true if there was a carry.
    #[must_use]
    fn slice_add_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool;

    /// Subtracts a value from the big integer slice, returning true if there was a borrow.
    #[must_use]
    fn slice_sub_value_assign(&mut self, value: Self::ValueT) -> bool;

    /// Multiplies the big integer slice by a value, returning any carry that results.
    #[must_use]
    fn slice_mul_value_assign(&mut self, value: Self::ValueT) -> Self::ValueT;

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

    /// Compares this big integer slice with another, returning an Ordering.
    #[must_use]
    fn slice_cmp(&self, other: &Self) -> Ordering;

    /// Adds another big integer slice to this one modulo a given modulus.
    fn slice_add_modulo_assign(&mut self, other: &Self, modulus: &Self);
}

impl<T: UnsignedInteger> BigIntegerOps for [T] {
    #[inline]
    fn slice_value_bits_count(&self) -> u32 {
        let mut bits_count = 0;
        for (i, &x) in self.iter().enumerate().rev() {
            if x.is_zero() {
                continue;
            }
            bits_count = (i as u32 + 1) * T::BITS - (x.leading_zeros());
            break;
        }

        bits_count
    }

    #[inline]
    fn slice_left_shift_assign(&mut self, bits: u32) {
        let mut pre = T::ZERO;
        let mut temp = T::ZERO;
        let right_shift_bits = T::BITS - bits;
        self.iter_mut().for_each(|value| {
            temp = *value;
            *value = *value << bits | pre >> right_shift_bits;
            pre = temp;
        });
    }

    #[inline]
    fn slice_right_shift_assign(&mut self, bits: u32) {
        let mut pre = T::ZERO;
        let mut temp = T::ZERO;
        let left_shift_bits = T::BITS - bits;
        self.iter_mut().rev().for_each(|value| {
            temp = *value;
            *value = pre << left_shift_bits | *value >> bits;
            pre = temp;
        });
    }

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

    fn slice_add_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool {
        let mut carry;

        let mut a_iter = self.iter();
        let mut b_iter = result.iter_mut();

        let a_first = a_iter.next().unwrap();
        let b_first = b_iter.next().unwrap();

        (*b_first, carry) = a_first.overflowing_add(value);

        for (a, b) in a_iter.zip(b_iter) {
            if !carry {
                return false;
            }
            (*b, carry) = a.overflowing_add(T::ONE);
        }

        carry
    }

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

    fn slice_add_assign(&mut self, other: &Self) -> bool {
        debug_assert_eq!(self.len(), other.len());

        let mut carry = false;

        for (xs, ys) in self.iter_mut().zip(other) {
            (*xs, carry) = xs.carrying_add(*ys, carry);
        }

        carry
    }

    fn slice_sub_assign(&mut self, other: &Self) -> bool {
        debug_assert_eq!(self.len(), other.len());

        let mut borrow = false;

        for (xs, ys) in self.iter_mut().zip(other) {
            (*xs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        borrow
    }

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

    fn slice_add_modulo_assign(&mut self, other: &Self, modulus: &Self) {
        debug_assert!(self.len() == other.len() && self.len() == modulus.len());
        debug_assert!(self.slice_cmp(modulus).is_lt());
        debug_assert!(other.slice_cmp(modulus).is_lt());

        let carry = self.slice_add_assign(other);
        if carry || self.slice_cmp(modulus).is_ge() {
            let _ = self.slice_sub_assign(modulus);
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
