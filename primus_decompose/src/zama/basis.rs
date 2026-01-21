use num_traits::ConstOne;
use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::{
    primitive::{ScalarIter, ValueMask},
    zama::common::SignedDecomposeIter,
};

/// The basis for approximate signed decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct ApproxSignedBasis<T: UnsignedInteger> {
    modulus: T,
    modulus_is_power_of_2: bool,
    basis: T,
    decompose_length: usize,
    value_bits: u32,
    decomposed_bits_minus_one: u32,
    log_basis: u32,
    drop_length: usize,
    scalars: Vec<T>,
    drop_value_masks: Vec<ValueMask<T>>,
    value_masks: Vec<ValueMask<T>>,
}

impl<T: UnsignedInteger> Eq for ApproxSignedBasis<T> {}

impl<T: UnsignedInteger> PartialEq for ApproxSignedBasis<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus
            && self.basis == other.basis
            && self.decompose_length == other.decompose_length
    }
}

impl<T: UnsignedInteger> ApproxSignedBasis<T> {
    #[inline]
    pub fn new(modulus: T, log_basis: u32, reverse_length: Option<usize>) -> Self {
        assert!(log_basis > 0 && T::BITS >= log_basis);

        let basis = <T as ConstOne>::ONE << log_basis;
        let basis_minus_one = basis - <T as ConstOne>::ONE;

        let modulus_is_power_of_2;
        let value_bits;

        if modulus.is_power_of_two() {
            modulus_is_power_of_2 = true;
            value_bits = modulus.trailing_zeros();
        } else {
            modulus_is_power_of_2 = false;
            value_bits = T::BITS - modulus.leading_zeros();
        }
        assert!(value_bits >= log_basis);

        let mut decompose_length = value_bits.div_ceil(log_basis) as usize;

        let mut drop_length = 0;

        if let Some(reverse_len) = reverse_length {
            assert!(decompose_length >= reverse_len);
            drop_length = decompose_length - reverse_len;
            decompose_length = reverse_len;
        }

        assert!(decompose_length > 0);

        let decomposed_bits_minus_one = (decompose_length + drop_length) as u32 * log_basis - 1;

        let mut scalars = vec![T::ZERO; decompose_length];
        let mut prev: Option<T> = None;
        scalars.iter_mut().for_each(|scalar| {
            if let Some(pre) = prev.as_mut() {
                *pre <<= log_basis;
                *scalar = *pre;
            } else {
                *scalar = T::ONE << (drop_length as u32 * log_basis);
                prev = Some(*scalar);
            }
        });

        let mut drop_value_masks = Vec::with_capacity(drop_length);
        let mut value_masks = Vec::with_capacity(decompose_length);
        let mut prev = ValueMask::new(basis_minus_one, 0);

        for _ in 0..drop_length {
            drop_value_masks.push(prev);
            prev = prev.next(log_basis);
        }

        for _ in 0..decompose_length {
            value_masks.push(prev);
            prev = prev.next(log_basis);
        }

        Self {
            modulus,
            modulus_is_power_of_2,
            basis,
            decompose_length,
            value_bits,
            decomposed_bits_minus_one,
            log_basis,
            drop_length,
            scalars,
            drop_value_masks,
            value_masks,
        }
    }

    /// Checks whether the modulus of this [`ApproxSignedBasis<T>`] is power of 2.
    #[inline]
    pub fn modulus_is_power_of_2(&self) -> bool {
        self.modulus_is_power_of_2
    }

    /// Returns the drop length of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn drop_length(&self) -> usize {
        self.drop_length
    }

    /// Returns the decompose length of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_length(&self) -> usize {
        self.decompose_length
    }

    /// Returns the basis value of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn basis_value(&self) -> T {
        self.basis
    }

    /// Returns the log basis of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn log_basis(&self) -> u32 {
        self.log_basis
    }

    /// Returns an iterator over the dropped signed decomposition operators of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn drop_iter<'a>(&'a self) -> SignedDecomposeIter<'a, T> {
        SignedDecomposeIter {
            value_masks: self.drop_value_masks.iter(),
            log_basis: self.log_basis,
            basis: self.basis,
            modulus: self.modulus,
        }
    }

    /// Returns an iterator over the signed decomposition operators of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_iter<'a>(&'a self) -> SignedDecomposeIter<'a, T> {
        SignedDecomposeIter {
            value_masks: self.value_masks.iter(),
            log_basis: self.log_basis,
            basis: self.basis,
            modulus: self.modulus,
        }
    }

    /// Returns an iterator over scalars of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn scalar_iter<'a>(&'a self) -> ScalarIter<'a, T> {
        ScalarIter::new(&self.scalars)
    }

    /// Init carry and adjusted value for a value.
    #[inline]
    pub fn init_value_sign(&self, value: T) -> (T, bool) {
        if value >> self.decomposed_bits_minus_one != T::ZERO {
            (self.modulus - value, true)
        } else {
            (value, false)
        }
    }

    pub fn decomposed_bits(&self) -> u32 {
        self.decomposed_bits_minus_one + 1
    }
}
