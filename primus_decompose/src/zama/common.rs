use std::iter::FusedIterator;

use primus_integer::UnsignedInteger;

use crate::primitive::ValueMask;

/// An iterator over the signed decomposition operators.
pub struct SignedDecomposeIter<'a, T: UnsignedInteger> {
    pub(super) value_masks: std::slice::Iter<'a, ValueMask<T>>,
    pub(super) log_basis: u32,
    pub(super) basis: T,
    pub(super) modulus: T,
}

impl<'a, T: UnsignedInteger> Iterator for SignedDecomposeIter<'a, T> {
    type Item = OnceSignedDecompose<T>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.value_masks
            .next()
            .map(|&value_mask| OnceSignedDecompose::<T> {
                value_mask,
                log_basis: self.log_basis,
                basis: self.basis,
                modulus: self.modulus,
            })
    }
}

impl<'a, T: UnsignedInteger> FusedIterator for SignedDecomposeIter<'a, T> {}

/// The signed decomposition operator which can execute once decomposition.
pub struct OnceSignedDecompose<T: UnsignedInteger> {
    value_mask: ValueMask<T>,
    log_basis: u32,
    basis: T,
    modulus: T,
}

impl<T: UnsignedInteger> OnceSignedDecompose<T> {
    /// Execute once decomposition and return the decomposed value and carry for next decomposition.
    #[inline]
    pub fn decompose(&self, value: T, sign: bool, carry: bool) -> (T, bool) {
        let mut temp = self.value_mask.get_value(value) + T::as_from(carry);

        let next_carry = !(temp >> self.log_basis).is_zero();
        if next_carry {
            temp -= self.basis
        }

        if sign && temp != T::ZERO {
            temp = self.modulus - temp;
        }

        (temp, next_carry)
    }

    /// Execute once decomposition, store carry for next decomposition back to `carry`.
    #[inline]
    pub fn decompose_inplace(
        &self,
        value: T,
        sign: bool,
        decomposed_value: &mut T,
        carry: &mut bool,
    ) {
        let mut temp = self.value_mask.get_value(value) + T::as_from(*carry);
        *carry = !(temp >> self.log_basis).is_zero();

        if *carry {
            temp -= self.basis
        }

        if sign && temp != T::ZERO {
            temp = self.modulus - temp;
        }

        *decomposed_value = temp;
    }
}
