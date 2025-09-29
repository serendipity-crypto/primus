use std::iter::FusedIterator;

use primus_integer::{BigIntegerOps, UnsignedInteger, izip};
use serde::{Deserialize, Serialize};

/// An iterator over scalars.
pub struct BigUintScalarIter<T: UnsignedInteger> {
    scalar: Vec<T>,
    length: usize,
    log_basis: u32,
}

impl<T: UnsignedInteger> BigUintScalarIter<T> {
    /// Creates a new [`BigUintScalarIter<T>`].
    #[inline]
    pub fn new(scalar: Vec<T>, length: usize, log_basis: u32) -> Self {
        Self {
            scalar,
            length,
            log_basis,
        }
    }
}

impl<T: UnsignedInteger> Iterator for BigUintScalarIter<T> {
    type Item = Vec<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.length == 0 {
            None
        } else {
            let next = self.scalar.clone();
            self.length -= 1;
            if self.length != 0 {
                self.scalar.slice_left_shift_assign(self.log_basis);
            }
            Some(next)
        }
    }
}

impl<T: UnsignedInteger> FusedIterator for BigUintScalarIter<T> {}

/// An iterator over the signed decomposition operators.
pub struct BigUintSignedDecomposeIter<T: UnsignedInteger> {
    pub(super) length: usize,
    pub(super) value_mask: ValueMask<T>,
    pub(super) mask_shl_bits: u32,
    pub(super) carry_mask: T,
    pub(super) basis_minus_one: T,
    pub(super) modulus_minus_basis: Vec<T>,
}

impl<T: UnsignedInteger> Iterator for BigUintSignedDecomposeIter<T> {
    type Item = OnceBigUintSignedDecompose<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.length == 0 {
            None
        } else {
            let next = OnceBigUintSignedDecompose::<T> {
                value_mask: self.value_mask,
                carry_mask: self.carry_mask,
                basis_minus_one: self.basis_minus_one,
                modulus_minus_basis: self.modulus_minus_basis.clone(),
            };

            self.length -= 1;

            if self.length > 0 {
                self.value_mask = self
                    .value_mask
                    .next(self.mask_shl_bits, self.basis_minus_one);
            }

            Some(next)
        }
    }
}

impl<T: UnsignedInteger> FusedIterator for BigUintSignedDecomposeIter<T> {}

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
        let mut mask = [basis_minus_one, T::ZERO];
        mask.slice_left_shift_assign(drop_bits);
        let shr_bits = drop_bits;
        let shl_bits = T::BITS - drop_bits;

        Self {
            index,
            mask,
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

        let mut mask = [basis_minus_one, T::ZERO];
        mask.slice_left_shift_assign(new_shr_bits);

        self.mask = mask;
        self.shr_bits = new_shr_bits;
        self.shl_bits = T::BITS - new_shr_bits;

        self
    }

    #[inline]
    fn get_value(&self, value: &[T]) -> T {
        let temp = (unsafe { *value.get_unchecked(self.index) } & self.mask[0]) >> self.shr_bits;

        if self.mask[1].is_zero() {
            temp
        } else {
            temp | (unsafe { *value.get_unchecked(self.index + 1) } & self.mask[1]) << self.shl_bits
        }
    }
}

/// The signed decomposition operator which can execute once decomposition.
pub struct OnceBigUintSignedDecompose<T: UnsignedInteger> {
    value_mask: ValueMask<T>,
    carry_mask: T,
    basis_minus_one: T,
    modulus_minus_basis: Vec<T>,
}

impl<T: UnsignedInteger> OnceBigUintSignedDecompose<T> {
    /// Execute once decomposition and return the decomposed value and carry for next decomposition.
    #[inline]
    pub fn decompose(&self, value: &[T], carry: bool) -> (Vec<T>, bool) {
        let temp = self.value_mask.get_value(value) + T::as_from(carry);

        let next_carry = !(temp & self.carry_mask).is_zero();
        let mut result = vec![T::ZERO; value.len()];
        if next_carry {
            if temp <= self.basis_minus_one {
                let _ = self
                    .modulus_minus_basis
                    .slice_add_value_inplace(temp, &mut result);
            }
        } else {
            result[0] = temp;
        }

        (result, next_carry)
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
                    .slice_add_value_inplace(temp, decomposed_value);
            }
        } else {
            decomposed_value.fill(T::ZERO);
            decomposed_value[0] = temp;
        }
    }

    /// Execute once decomposition for slice, store carries for next decomposition back to `carries`.
    #[inline]
    pub fn decompose_slice_inplace(
        &self,
        values: &[T],
        decomposed_values: &mut [T],
        carries: &mut [bool],
        value_chunk_size: usize,
    ) {
        for (value, decomposed_value, carry) in izip!(
            values.chunks_exact(value_chunk_size),
            decomposed_values.chunks_exact_mut(value_chunk_size),
            carries
        ) {
            self.decompose_inplace(value, decomposed_value, carry);
        }
    }
}
