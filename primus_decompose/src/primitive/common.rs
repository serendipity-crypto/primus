use std::iter::FusedIterator;

use primus_integer::{UnsignedInteger, izip};
use serde::{Deserialize, Serialize};

/// An iterator over scalars.
pub struct ScalarIter<'a, T: UnsignedInteger> {
    iter: std::iter::Copied<std::slice::Iter<'a, T>>,
}

impl<'a, T: UnsignedInteger> ScalarIter<'a, T> {
    /// Creates a new [`ScalarIter<T>`].
    #[inline]
    pub fn new(scalars: &'a [T]) -> Self {
        Self {
            iter: scalars.iter().copied(),
        }
    }
}

impl<'a, T: UnsignedInteger> Iterator for ScalarIter<'a, T> {
    type Item = T;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<'a, T: UnsignedInteger> FusedIterator for ScalarIter<'a, T> {}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct ValueMask<T: UnsignedInteger> {
    mask: T,
    shr_bits: u32,
}

impl<T: UnsignedInteger> ValueMask<T> {
    #[inline]
    pub fn new(basis_minus_one: T, drop_bits: u32) -> Self {
        Self {
            mask: basis_minus_one << drop_bits,
            shr_bits: drop_bits,
        }
    }

    #[inline]
    pub fn next(self, mask_shl_bits: u32) -> Self {
        Self {
            mask: self.mask << mask_shl_bits,
            shr_bits: self.shr_bits + mask_shl_bits,
        }
    }

    #[inline]
    fn get_value(&self, value: T) -> T {
        (value & self.mask) >> self.shr_bits
    }
}

/// An iterator over the signed decomposition operators.
pub struct SignedDecomposeIter<'a, T: UnsignedInteger> {
    pub(super) value_masks: std::slice::Iter<'a, ValueMask<T>>,
    pub(super) carry_mask: T,
    pub(super) basis_minus_one: T,
    pub(super) modulus_minus_basis: T,
}

impl<'a, T: UnsignedInteger> Iterator for SignedDecomposeIter<'a, T> {
    type Item = OnceSignedDecompose<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.value_masks
            .next()
            .map(|&value_mask| OnceSignedDecompose::<T> {
                value_mask,
                carry_mask: self.carry_mask,
                basis_minus_one: self.basis_minus_one,
                modulus_minus_basis: self.modulus_minus_basis,
            })
    }
}

impl<'a, T: UnsignedInteger> FusedIterator for SignedDecomposeIter<'a, T> {}

/// The signed decomposition operator which can execute once decomposition.
pub struct OnceSignedDecompose<T: UnsignedInteger> {
    value_mask: ValueMask<T>,
    carry_mask: T,
    basis_minus_one: T,
    modulus_minus_basis: T,
}

impl<T: UnsignedInteger> OnceSignedDecompose<T> {
    /// Execute once decomposition and return the decomposed value and carry for next decomposition.
    #[inline]
    pub fn decompose(&self, value: T, carry: bool) -> (T, bool) {
        let mut temp = self.value_mask.get_value(value) + T::as_from(carry);

        let next_carry = !(temp & self.carry_mask).is_zero();
        if next_carry {
            if temp > self.basis_minus_one {
                temp = T::ZERO;
            } else {
                temp += self.modulus_minus_basis
            }
        }

        (temp, next_carry)
    }

    /// Execute once decomposition, store carry for next decomposition back to `carry`.
    #[inline]
    pub fn decompose_inplace(&self, value: T, decomposed_value: &mut T, carry: &mut bool) {
        let temp = self.value_mask.get_value(value) + T::as_from(*carry);
        *carry = !(temp & self.carry_mask).is_zero();

        if *carry {
            if temp > self.basis_minus_one {
                *decomposed_value = T::ZERO;
            } else {
                *decomposed_value = temp + self.modulus_minus_basis
            }
        } else {
            *decomposed_value = temp;
        }
    }

    /// Execute once decomposition for slice, store carries for next decomposition back to `carries`.
    #[inline]
    pub fn decompose_slice_inplace(
        &self,
        values: &[T],
        decomposed_values: &mut [T],
        carries: &mut [bool],
    ) {
        for (&value, decomposed_value, carry) in izip!(values, decomposed_values, carries) {
            self.decompose_inplace(value, decomposed_value, carry);
        }
    }
}
