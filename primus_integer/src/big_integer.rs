use core::cmp::Ordering;

use crate::{Data, DataMut, RawData, UnsignedInteger, impl_iters, izip};

#[repr(transparent)]
pub struct BigUint<S>(pub S)
where
    S: RawData,
    <S as RawData>::Elem: UnsignedInteger;

impl_iters!(BigUint, bit_uint);

impl<S> Clone for BigUint<S>
where
    S: RawData + Clone,
    <S as RawData>::Elem: UnsignedInteger,
{
    #[inline(always)]
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<S> Copy for BigUint<S>
where
    S: RawData + Copy,
    <S as RawData>::Elem: UnsignedInteger,
{
}

impl<S, A, T> PartialEq<BigUint<A>> for BigUint<S>
where
    S: Data<Elem = T>,
    A: Data<Elem = T>,
    T: UnsignedInteger,
{
    #[inline]
    fn eq(&self, other: &BigUint<A>) -> bool {
        assert_eq!(self.len(), other.len());
        self.iter().zip(other.iter()).all(|(&a, &b)| a == b)
    }
}

impl<S, T> BigUint<S>
where
    S: Data<Elem = T>,
    T: UnsignedInteger,
{
    #[allow(clippy::len_without_is_empty)]
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline(always)]
    pub fn digits(&self) -> &[T] {
        self.0.as_slice()
    }

    #[inline(always)]
    pub fn iter<'a>(&'a self) -> std::slice::Iter<'a, T> {
        self.0.iter()
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.iter().all(T::is_zero)
    }

    /// Gets the bits count of the big unsigned integer.
    #[must_use]
    #[inline]
    pub fn bits_count(&self) -> u32 {
        self.iter()
            .enumerate()
            .rev()
            .find(|(_, v)| !v.is_zero())
            .map_or(0, |(i, v)| T::BITS * (i as u32 + 1) - v.leading_zeros())
    }

    /// Adds a value to the big integer, returning true if there was a carry.
    #[must_use]
    #[inline]
    pub fn add_value_inplace<A>(&self, value: T, result: &mut BigUint<A>) -> bool
    where
        A: DataMut<Elem = T>,
    {
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

    /// Subtracts a value to the big integer, returning true if there was a borrow.
    #[must_use]
    #[inline]
    pub fn sub_value_inplace<A>(&self, value: T, result: &mut BigUint<A>) -> bool
    where
        A: DataMut<Elem = T>,
    {
        debug_assert_eq!(self.len(), result.len());

        let mut borrow;

        let mut a_iter = self.iter();
        let mut b_iter = result.iter_mut();

        let a_first = a_iter.next().unwrap();
        let b_first = b_iter.next().unwrap();

        (*b_first, borrow) = a_first.overflowing_sub(value);

        while borrow {
            if let Some(a_next) = a_iter.next()
                && let Some(b_next) = b_iter.next()
            {
                (*b_next, borrow) = a_next.overflowing_sub(T::ONE);
            } else {
                return borrow;
            }
        }

        for (a, b) in a_iter.zip(b_iter) {
            *b = *a;
        }

        borrow
    }

    /// Multiplies the big integer by a value, storing the result in another big integer.
    #[must_use]
    #[inline]
    pub fn mul_value_inplace<A>(&self, value: T, result: &mut BigUint<A>) -> T
    where
        A: DataMut<Elem = T>,
    {
        debug_assert_eq!(result.len(), self.len());

        if value.is_zero() {
            result.set_zero();
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for (ele, res) in self.iter().zip(result.iter_mut()) {
            (*res, carry) = value.carrying_mul(*ele, carry);
        }

        carry
    }

    /// Multiplies the big integer by a value, then add to another big integer.
    #[must_use]
    #[inline]
    pub fn mul_value_add_inplace<A>(&self, value: T, result: &mut BigUint<A>) -> T
    where
        A: DataMut<Elem = T>,
    {
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

    /// Adds two big integers to the result, returning true if there was a carry.
    #[must_use]
    #[inline]
    pub fn add_inplace<A, B>(&self, other: &BigUint<A>, result: &mut BigUint<B>) -> bool
    where
        A: Data<Elem = T>,
        B: DataMut<Elem = T>,
    {
        debug_assert_eq!(self.len(), other.len());
        debug_assert_eq!(self.len(), result.len());

        let mut carry = false;
        for (xs, ys, zs) in izip!(self.iter(), other.iter(), result.iter_mut()) {
            (*zs, carry) = xs.carrying_add(*ys, carry);
        }

        carry
    }

    /// Subtracts another big integer from this one, returning true if there was a borrow.
    #[must_use]
    #[inline]
    pub fn sub_inplace<A, B>(&self, other: &BigUint<A>, result: &mut BigUint<B>) -> bool
    where
        A: Data<Elem = T>,
        B: DataMut<Elem = T>,
    {
        debug_assert_eq!(self.len(), other.len());
        debug_assert_eq!(self.len(), result.len());

        let mut borrow = false;
        for (xs, ys, zs) in izip!(self.iter(), other.iter(), result.iter_mut()) {
            (*zs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        borrow
    }

    /// Compares this big integer with another, returning an [`Ordering`].
    #[must_use]
    #[inline]
    pub fn cmp<A>(&self, other: &BigUint<A>) -> Ordering
    where
        A: Data<Elem = T>,
    {
        debug_assert_eq!(self.len(), other.len());

        for (a, b) in self.iter().rev().zip(other.iter().rev()) {
            match a.cmp(b) {
                Ordering::Equal => continue,
                neq => return neq,
            }
        }

        Ordering::Equal
    }

    /// Adds two big integers to the result modulo a given modulus.
    #[inline]
    pub fn add_modulo_inplace<A, B, C>(
        &self,
        other: &BigUint<A>,
        result: &mut BigUint<B>,
        modulus: &BigUint<C>,
    ) where
        A: Data<Elem = T>,
        B: DataMut<Elem = T>,
        C: Data<Elem = T>,
    {
        debug_assert!(
            self.len() == other.len() && self.len() == result.len() && self.len() == modulus.len()
        );
        debug_assert!(self.cmp(modulus).is_lt());
        debug_assert!(other.cmp(modulus).is_lt());

        let carry = self.add_inplace(other, result);
        if carry || result.cmp(modulus).is_ge() {
            let _ = result.sub_assign(modulus);
        }
    }

    /// Subs another big integer to this one modulo a given modulus.
    #[inline]
    pub fn sub_modulo_inplace<A, B, C>(
        &self,
        other: &BigUint<A>,
        result: &mut BigUint<B>,
        modulus: &BigUint<C>,
    ) where
        A: Data<Elem = T>,
        B: DataMut<Elem = T>,
        C: Data<Elem = T>,
    {
        debug_assert!(
            self.len() == other.len() && self.len() == result.len() && self.len() == modulus.len()
        );
        debug_assert!(self.cmp(modulus).is_lt());
        debug_assert!(other.cmp(modulus).is_lt());

        if self.sub_inplace(other, result) {
            let _ = result.add_assign(modulus);
        }
    }

    /// Negates the big integer modulo a given modulus.
    #[inline]
    pub fn neg_modulo_inplace<A, B>(&self, result: &mut BigUint<A>, modulus: &BigUint<B>)
    where
        A: DataMut<Elem = T>,
        B: Data<Elem = T>,
    {
        debug_assert!(self.len() == result.len() && self.len() == modulus.len());
        debug_assert!(self.cmp(modulus).is_lt());

        if self.is_zero() {
            result.set_zero();
        } else {
            let mut borrow = false;
            for (xs, ys, zs) in izip!(self.iter(), modulus.iter(), result.iter_mut()) {
                (*zs, borrow) = ys.borrowing_sub(*xs, borrow);
            }
        }
    }
}

impl<S, T> BigUint<S>
where
    S: DataMut<Elem = T>,
    T: UnsignedInteger,
{
    #[inline(always)]
    pub fn digits_mut(&mut self) -> &mut [T] {
        self.0.as_mut_slice()
    }

    #[inline(always)]
    pub fn iter_mut<'a>(&'a mut self) -> std::slice::IterMut<'a, T> {
        self.0.iter_mut()
    }

    #[inline(always)]
    pub fn set_zero(&mut self) {
        self.0.fill(T::ZERO);
    }

    /// Left shifts the big integer.
    #[inline]
    pub fn left_shift_assign(&mut self, bits: u32) -> T {
        if bits != 0 {
            let mut pre = T::ZERO;
            let mut temp = T::ZERO;
            let right_shift_bits = T::BITS - bits;
            self.iter_mut().for_each(|value| {
                temp = *value;
                *value = *value << bits | pre >> right_shift_bits;
                pre = temp;
            });
            pre >> right_shift_bits
        } else {
            T::ZERO
        }
    }

    /// Right shifts the big integer.
    #[inline]
    pub fn right_shift_assign(&mut self, bits: u32) {
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

    /// Adds a value to the big integer, returning true if there was a carry.
    #[must_use]
    #[inline]
    pub fn add_value_assign(&mut self, value: T) -> bool {
        let mut carry;
        match self.digits_mut() {
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

    /// Subtracts a value from the big integer, returning true if there was a borrow.
    #[must_use]
    #[inline]
    pub fn sub_value_assign(&mut self, value: T) -> bool {
        let mut borrow;
        match self.digits_mut() {
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

    /// Multiplies the big integer by a value, returning any carry that results.
    #[must_use]
    #[inline]
    pub fn mul_value_assign(&mut self, value: T) -> T {
        if value.is_zero() {
            self.set_zero();
            return T::ZERO;
        }

        let mut carry = T::ZERO;
        for ele in self.iter_mut() {
            (*ele, carry) = value.carrying_mul(*ele, carry);
        }

        carry
    }

    /// Adds another big integer to this one, returning true if there was a carry.
    #[must_use]
    #[inline]
    pub fn add_assign<A>(&mut self, other: &BigUint<A>) -> bool
    where
        A: Data<Elem = T>,
    {
        debug_assert_eq!(self.len(), other.len());

        let mut carry = false;

        for (xs, ys) in self.iter_mut().zip(other.iter()) {
            (*xs, carry) = xs.carrying_add(*ys, carry);
        }

        carry
    }

    /// Subtracts another big integer from this one, returning true if there was a borrow.
    #[must_use]
    #[inline]
    pub fn sub_assign<A>(&mut self, other: &BigUint<A>) -> bool
    where
        A: Data<Elem = T>,
    {
        debug_assert_eq!(self.len(), other.len());

        let mut borrow = false;

        for (xs, ys) in self.iter_mut().zip(other.iter()) {
            (*xs, borrow) = xs.borrowing_sub(*ys, borrow);
        }

        borrow
    }

    /// Adds another big integer to this one modulo a given modulus.
    #[inline]
    pub fn add_modulo_assign<A, B>(&mut self, other: &BigUint<A>, modulus: &BigUint<B>)
    where
        A: Data<Elem = T>,
        B: Data<Elem = T>,
    {
        debug_assert!(self.len() == other.len() && self.len() == modulus.len());
        debug_assert!(self.cmp(modulus).is_lt());
        debug_assert!(other.cmp(modulus).is_lt());

        let carry = self.add_assign(other);
        if carry || self.cmp(modulus).is_ge() {
            let _ = self.sub_assign(modulus);
        }
    }

    /// Subs another big integer from this one modulo a given modulus.
    #[inline]
    pub fn sub_modulo_assign<A, B>(&mut self, other: &BigUint<A>, modulus: &BigUint<B>)
    where
        A: Data<Elem = T>,
        B: Data<Elem = T>,
    {
        debug_assert!(self.len() == other.len() && self.len() == modulus.len());
        debug_assert!(self.cmp(modulus).is_lt());
        debug_assert!(other.cmp(modulus).is_lt());

        if self.sub_assign(other) {
            let _ = self.add_assign(modulus);
        }
    }

    /// Negates the big integer modulo a given modulus.
    #[inline]
    pub fn neg_modulo_assign<A>(&mut self, modulus: &BigUint<A>)
    where
        A: Data<Elem = T>,
    {
        debug_assert!(self.len() == modulus.len());
        debug_assert!(self.cmp(modulus).is_lt());

        if !self.is_zero() {
            let mut borrow = false;
            for (xs, ys) in self.iter_mut().zip(modulus.iter()) {
                (*xs, borrow) = ys.borrowing_sub(*xs, borrow);
            }
        }
    }
}

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
    #[must_use]
    fn slice_left_shift_assign(&mut self, bits: u32) -> Self::ValueT;

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

    /// Subtracts a value to the big integer slice, returning true if there was a borrow.
    #[must_use]
    fn slice_sub_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool;

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

    /// Negates the big integer slice modulo a given modulus.
    fn slice_neg_modulo_assign(&mut self, modulus: &Self);

    /// Adds two big integer slices to result modulo a given modulus.
    fn slice_add_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self);

    /// Subs another big integer slice to this one modulo a given modulus.
    fn slice_sub_modulo_inplace(&self, other: &Self, result: &mut Self, modulus: &Self);

    /// Negates the big integer slice modulo a given modulus.
    fn slice_neg_modulo_inplace(&self, result: &mut Self, modulus: &Self);
}

impl<T: UnsignedInteger> BigIntegerOps for [T] {
    #[inline]
    fn slice_left_shift_assign(&mut self, bits: u32) -> Self::ValueT {
        if bits != 0 {
            let mut pre = T::ZERO;
            let mut temp = T::ZERO;
            let right_shift_bits = T::BITS - bits;
            self.iter_mut().for_each(|value| {
                temp = *value;
                *value = *value << bits | pre >> right_shift_bits;
                pre = temp;
            });
            pre >> right_shift_bits
        } else {
            T::ZERO
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

    fn slice_sub_value_inplace(&self, value: Self::ValueT, result: &mut Self) -> bool {
        debug_assert_eq!(self.len(), result.len());

        let mut borrow;

        let mut a_iter = self.iter();
        let mut b_iter = result.iter_mut();

        let a_first = a_iter.next().unwrap();
        let b_first = b_iter.next().unwrap();

        (*b_first, borrow) = a_first.overflowing_sub(value);

        while borrow {
            if let Some(a_next) = a_iter.next()
                && let Some(b_next) = b_iter.next()
            {
                (*b_next, borrow) = a_next.overflowing_sub(T::ONE);
            } else {
                return borrow;
            }
        }

        for (a, b) in a_iter.zip(b_iter) {
            *b = *a;
        }

        borrow
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
    let mut result = BigUint(Vec::with_capacity(values.len()));
    result.0.push(values[0]);
    for &v in values.iter().skip(1) {
        let carry = result.mul_value_assign(v);
        if !carry.is_zero() {
            result.0.push(carry);
        }
    }
    result.0.shrink_to_fit();
    result.0
}

/// Multiplies many values together, except for one specified by index, returning the result as a big integer slice.
pub fn multiply_many_values_except<T: UnsignedInteger>(values: &[T], except: usize) -> Vec<T> {
    let mut result = BigUint(Vec::with_capacity(values.len() - 1));
    result.0.push(T::ONE);

    for (i, &v) in values.iter().enumerate() {
        if i == except {
            continue;
        }
        let carry = result.mul_value_assign(v);
        if !carry.is_zero() {
            result.0.push(carry);
        }
    }

    result.0.shrink_to_fit();
    result.0
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
        let carry = BigUint(&mut result[0..len]).mul_value_assign(v);
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
        for &r in value.iter().rev() {
            result <<= ValueT::BITS;
            result |= r as u128;
        }
        result
    }

    #[test]
    fn test_big_uint_ops() {
        let mut rng = rand::rng();
        let moduli: [ValueT; 3] = [134215681, 134176769, 132120577];
        let modulus = BigUint(multiply_many_values(&moduli));
        let m_raw = compose(modulus.digits());

        assert_eq!(128 - m_raw.leading_zeros(), modulus.bits_count());

        let distrs = moduli.map(|m| Uniform::new(0, m).unwrap());

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a = BigUint(multiply_many_values(&a_residues));
        let mut a_raw = compose(a.digits());

        a.right_shift_assign(3);
        a_raw >>= 3;
        assert_eq!(a_raw, compose(a.digits()));

        let carry = a.left_shift_assign(3);
        assert_eq!(carry, 0);
        a_raw <<= 3;
        assert_eq!(a_raw, compose(a.digits()));

        let v: ValueT = rng.random();
        let _r = a.add_value_assign(v);
        a_raw += v as u128;
        assert_eq!(a_raw, compose(a.digits()));

        let _r = a.sub_value_assign(v);
        a_raw -= v as u128;
        assert_eq!(a_raw, compose(a.digits()));

        let r = a.mul_value_assign(v);
        let mut p = a.clone();
        p.0.push(r);
        a_raw *= v as u128;
        assert_eq!(a_raw, compose(p.digits()));

        let mut result = BigUint(vec![0; a.len()]);
        a_raw = compose(a.digits());
        let _carry = a.add_value_inplace(v, &mut result);
        assert_eq!(a_raw + v as u128, compose(result.digits()));

        let _borrow = a.sub_value_inplace(v, &mut result);
        assert_eq!(a_raw - v as u128, compose(result.digits()));

        let r = a.mul_value_inplace(v, &mut result);
        result.0.push(r);
        assert_eq!(a_raw * v as u128, compose(result.digits()));

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let b_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a = BigUint(multiply_many_values(&a_residues));
        let b = BigUint(multiply_many_values(&b_residues));
        let a_raw = compose(a.digits());
        let b_raw = compose(b.digits());

        let mut result = b.clone();
        let carry = a.mul_value_add_inplace(v, &mut result);
        result.0.push(carry);
        assert_eq!(a_raw * v as u128 + b_raw, compose(result.digits()));

        let _r = a.add_assign(&b);
        assert_eq!(a_raw + b_raw, compose(a.digits()));

        let _r = a.sub_assign(&b);
        assert_eq!(a_raw, compose(a.digits()));

        a.add_modulo_assign(&b, &modulus);
        let r = (a_raw + b_raw) % m_raw;
        assert_eq!(r, compose(a.digits()));

        let a_residues = distrs.map(|distr| distr.sample(&mut rng));
        let b_residues = distrs.map(|distr| distr.sample(&mut rng));
        let mut a = BigUint(multiply_many_values(&a_residues));
        let b = BigUint(multiply_many_values(&b_residues));
        let a_raw = compose(a.digits());
        let b_raw = compose(b.digits());

        a.sub_modulo_assign(&b, &modulus);
        let r = (a_raw + m_raw - b_raw) % m_raw;
        assert_eq!(r, compose(a.digits()))
    }
}
