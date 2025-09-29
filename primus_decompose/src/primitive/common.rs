use std::iter::FusedIterator;

use primus_integer::{UnsignedInteger, izip};

/// An iterator over scalars.
pub struct ScalarIter<T: UnsignedInteger> {
    scalar: T,
    length: usize,
    log_basis: u32,
}

impl<T: UnsignedInteger> ScalarIter<T> {
    /// Creates a new [`ScalarIter<T>`].
    #[inline]
    pub fn new(scalar: T, length: usize, log_basis: u32) -> Self {
        Self {
            scalar,
            length,
            log_basis,
        }
    }
}

impl<T: UnsignedInteger> Iterator for ScalarIter<T> {
    type Item = T;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.length == 0 {
            None
        } else {
            let next = self.scalar;
            self.length -= 1;
            if self.length != 0 {
                self.scalar <<= self.log_basis;
            }
            Some(next)
        }
    }
}

impl<T: UnsignedInteger> FusedIterator for ScalarIter<T> {}

/// An iterator over the signed decomposition operators.
pub struct SignedDecomposeIter<T: UnsignedInteger> {
    pub(super) length: usize,
    pub(super) value_chunk_mask: T,
    pub(super) mask_shl_bits: u32,
    pub(super) value_shr_bits: u32,
    pub(super) carry_mask: T,
    pub(super) basis_minus_one: T,
    pub(super) modulus_minus_basis: T,
}

impl<T: UnsignedInteger> Iterator for SignedDecomposeIter<T> {
    type Item = OnceSignedDecompose<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.length == 0 {
            None
        } else {
            let next = OnceSignedDecompose::<T> {
                value_chunk_mask: self.value_chunk_mask,
                shr_bits: self.value_shr_bits,
                carry_mask: self.carry_mask,
                basis_minus_one: self.basis_minus_one,
                modulus_minus_basis: self.modulus_minus_basis,
            };

            self.length -= 1;

            if self.length > 0 {
                self.value_chunk_mask <<= self.mask_shl_bits;
                self.value_shr_bits += self.mask_shl_bits;
            }

            Some(next)
        }
    }
}

impl<T: UnsignedInteger> FusedIterator for SignedDecomposeIter<T> {}

/// The signed decomposition operator which can execute once decomposition.
pub struct OnceSignedDecompose<T: UnsignedInteger> {
    value_chunk_mask: T,
    shr_bits: u32,
    carry_mask: T,
    basis_minus_one: T,
    modulus_minus_basis: T,
}

impl<T: UnsignedInteger> OnceSignedDecompose<T> {
    /// Execute once decomposition and return the decomposed value and carry for next decomposition.
    #[inline]
    pub fn decompose(&self, value: T, carry: bool) -> (T, bool) {
        let mut temp = ((value & self.value_chunk_mask) >> self.shr_bits) + T::as_from(carry);

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
        let temp = ((value & self.value_chunk_mask) >> self.shr_bits) + T::as_from(*carry);
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
