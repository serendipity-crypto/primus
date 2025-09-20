use core::cmp::Ordering;

use integer::UnsignedInteger;

pub trait BigIntOps {
    type ValueT;
    fn slice_add_value_assign(&mut self, value: Self::ValueT) -> bool;
    fn slice_sub_value_assign(&mut self, value: Self::ValueT) -> bool;
    fn slice_mul_value_assign(&mut self, value: Self::ValueT) -> Self::ValueT;

    fn slice_mul_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT;

    fn slice_add_assign(&mut self, other: &Self) -> bool;
    fn slice_sub_assign(&mut self, other: &Self) -> bool;

    fn slice_cmp(&self, other: &Self) -> Ordering;

    fn slice_add_modulo_assign(&mut self, other: &Self, modulus: &Self);
}

impl<T: UnsignedInteger> BigIntOps for [T] {
    type ValueT = T;

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
            let (low, high) = value.carrying_mul(*ele, carry);
            *ele = low;
            carry = high;
        }

        carry
    }

    fn slice_mul_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> Self::ValueT {
        debug_assert!(result.len() >= self.len());

        if value.is_zero() {
            result.fill(T::ZERO);
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for (ele, res) in self.iter().zip(result.iter_mut()) {
            let (low, high) = value.carrying_mul(*ele, carry);
            *res = low;
            carry = high;
        }

        if result.len() > self.len() {
            result[self.len()..].fill(T::ZERO);
            result[self.len()] = carry;
            T::ZERO
        } else {
            carry
        }
    }

    fn slice_add_assign(&mut self, other: &Self) -> bool {
        debug_assert!(self.len() >= other.len());

        let mut carry = false;

        for (xs, ys) in self.iter_mut().zip(other) {
            (*xs, carry) = xs.carrying_add(*ys, carry);
        }

        if carry && self.len() > other.len() {
            for xs in self[other.len()..].iter_mut() {
                (*xs, carry) = xs.overflowing_add(T::ONE);
                if !carry {
                    return false;
                }
            }
        }

        carry
    }

    fn slice_sub_assign(&mut self, other: &Self) -> bool {
        debug_assert!(self.len() >= other.len());

        let mut borrow = false;

        for (xs, ys) in self[..other.len()].iter_mut().zip(other) {
            (*xs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        if borrow && self.len() > other.len() {
            for xs in self[other.len()..].iter_mut() {
                (*xs, borrow) = xs.overflowing_sub(T::ONE);
                if !borrow {
                    return false;
                }
            }
        }

        borrow
    }

    fn slice_cmp(&self, other: &Self) -> Ordering {
        assert_eq!(self.len(), other.len());

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
            self.slice_sub_assign(modulus);
        }
    }
}

pub fn multiply_many_values<T: UnsignedInteger>(values: &[T]) -> Vec<T> {
    let mut result = Vec::with_capacity(values.len());
    result.push(values[0]);
    for &v in values.iter().skip(1) {
        let carry = result.as_mut_slice().slice_mul_value_assign(v);
        if !carry.is_zero() {
            result.push(carry);
        }
    }
    result.shrink_to_fit();
    result
}

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
