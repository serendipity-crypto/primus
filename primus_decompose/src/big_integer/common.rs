use std::iter::FusedIterator;

use primus_integer::{BigUint, UnsignedInteger, izip};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct ValueMask<T: UnsignedInteger> {
    index: usize,
    mask: [T; 2],
    shr_bits: u32,
    shl_bits: u32,
}

impl<T: UnsignedInteger> ValueMask<T> {
    #[inline]
    pub fn new(basis_minus_one: T, drop_bits: u32) -> Self {
        let index = (drop_bits / T::BITS) as usize;
        let drop_bits = drop_bits % T::BITS;
        let mut mask = BigUint([basis_minus_one, T::ZERO]);
        let carry = mask.left_shift_assign(drop_bits);
        assert_eq!(carry, T::ZERO);
        let shr_bits = drop_bits;
        let shl_bits = T::BITS - drop_bits;

        Self {
            index,
            mask: mask.0,
            shr_bits,
            shl_bits,
        }
    }

    #[inline]
    pub fn next(mut self, mask_shl_bits: u32, basis_minus_one: T) -> Self {
        let mut new_shr_bits = mask_shl_bits + self.shr_bits;

        if new_shr_bits >= T::BITS {
            self.index += 1;
            new_shr_bits -= T::BITS;
        }

        let mut mask = BigUint([basis_minus_one, T::ZERO]);
        let carry = mask.left_shift_assign(new_shr_bits);
        assert_eq!(carry, T::ZERO);

        self.mask = mask.0;
        self.shr_bits = new_shr_bits;
        self.shl_bits = T::BITS - new_shr_bits;

        self
    }

    #[inline]
    fn get_value(&self, value: &[T]) -> T {
        let temp = (value[self.index] & self.mask[0]) >> self.shr_bits;

        if self.mask[1].is_zero() {
            temp
        } else {
            temp | (value[self.index + 1] & self.mask[1]) << self.shl_bits
        }
    }
}

pub struct BigUintSignedDecomposerIter<'a, T: UnsignedInteger> {
    pub(super) value_masks: std::slice::Iter<'a, ValueMask<T>>,
    pub(super) carry_mask: T,
    pub(super) basis_minus_one: T,
    pub(super) modulus_minus_basis: &'a [T],
}

impl<'a, T: UnsignedInteger> Iterator for BigUintSignedDecomposerIter<'a, T> {
    type Item = OnceBigUintSignedDecomposer<'a, T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.value_masks
            .next()
            .map(|&value_mask| OnceBigUintSignedDecomposer {
                value_mask,
                carry_mask: self.carry_mask,
                basis_minus_one: self.basis_minus_one,
                modulus_minus_basis: BigUint(self.modulus_minus_basis),
            })
    }
}

impl<'a, T: UnsignedInteger> FusedIterator for BigUintSignedDecomposerIter<'a, T> {}

/// The signed decomposition operator which can execute once decomposition.
pub struct OnceBigUintSignedDecomposer<'a, T: UnsignedInteger> {
    pub(super) value_mask: ValueMask<T>,
    pub(super) carry_mask: T,
    pub(super) basis_minus_one: T,
    pub(super) modulus_minus_basis: BigUint<&'a [T]>,
}

impl<'a, T: UnsignedInteger> OnceBigUintSignedDecomposer<'a, T> {
    /// Execute once decomposition and return the decomposed value and carry for next decomposition.
    #[inline]
    pub fn decompose(&self, value: &[T], carry: bool) -> (Vec<T>, bool) {
        let temp = self.value_mask.get_value(value) + T::as_from(carry);

        let next_carry = !(temp & self.carry_mask).is_zero();
        let mut result = BigUint(vec![T::ZERO; value.len()]);
        if next_carry {
            if temp <= self.basis_minus_one {
                let _ = self
                    .modulus_minus_basis
                    .add_value_inplace(temp, &mut result);
            }
        } else {
            result[0] = temp;
        }

        (result.0, next_carry)
    }

    #[inline]
    pub fn unsigned_decompose(&self, value: &[T], carry: bool) -> (T, bool) {
        let temp = self.value_mask.get_value(value) + T::as_from(carry);

        let next_carry = !(temp & self.carry_mask).is_zero();

        (temp & self.basis_minus_one, next_carry)
    }

    /// Execute once decomposition, store carry for next decomposition back to `carry`.
    #[inline]
    pub fn decompose_inplace(&self, value: &[T], decomposed_value: &mut [T], carry: &mut bool) {
        let temp = self.value_mask.get_value(value) + T::as_from(*carry);
        *carry = !(temp & self.carry_mask).is_zero();

        if *carry {
            if temp > self.basis_minus_one {
                decomposed_value.fill(T::ZERO);
            } else {
                let _ = self
                    .modulus_minus_basis
                    .add_value_inplace(temp, &mut BigUint(decomposed_value));
            }
        } else {
            decomposed_value.fill(T::ZERO);
            decomposed_value[0] = temp;
        }
    }

    /// Execute once decomposition, store carry for next decomposition back to `carry`.
    #[inline]
    pub fn unsigned_decompose_inplace(
        &self,
        value: &[T],
        decomposed_unsigned_value: &mut T,
        carry: &mut bool,
    ) {
        let temp = self.value_mask.get_value(value) + T::as_from(*carry);
        *carry = !(temp & self.carry_mask).is_zero();

        *decomposed_unsigned_value = temp & self.basis_minus_one;
    }

    /// Execute once decomposition for slice, store carries for next decomposition back to `carries`.
    #[inline]
    pub fn decompose_slice_inplace(
        &self,
        big_uint_values: &[T],
        decomposed_big_uint_values: &mut [T],
        carries: &mut [bool],
        big_uint_value_len: usize,
    ) {
        debug_assert_eq!(decomposed_big_uint_values.len(), big_uint_values.len());
        debug_assert_eq!(big_uint_values.len(), carries.len() * big_uint_value_len);
        for (value, decomposed_value, carry) in izip!(
            big_uint_values.chunks_exact(big_uint_value_len),
            decomposed_big_uint_values.chunks_exact_mut(big_uint_value_len),
            carries
        ) {
            self.decompose_inplace(value, decomposed_value, carry);
        }
    }

    /// Execute once decomposition for slice, store carries for next decomposition back to `carries`.
    #[inline]
    pub fn unsigned_decompose_slice_inplace(
        &self,
        big_uint_values: &[T],
        decomposed_unsigned_values: &mut [T],
        carries: &mut [bool],
        big_uint_value_len: usize,
    ) {
        debug_assert_eq!(carries.len(), decomposed_unsigned_values.len());
        debug_assert_eq!(big_uint_values.len(), carries.len() * big_uint_value_len);
        for (value, decomposed_unsigned_value, carry) in izip!(
            big_uint_values.chunks_exact(big_uint_value_len),
            decomposed_unsigned_values.iter_mut(),
            carries
        ) {
            self.unsigned_decompose_inplace(value, decomposed_unsigned_value, carry);
        }
    }
}
