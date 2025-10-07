use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::NttRlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different powers
/// of a base, used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct NttRlev<T: UnsignedInteger> {
    /// A vector of RLWE ciphertexts, each encrypted message with a different `basis`.
    data: Vec<NttRlwe<T>>,
}

impl<T: UnsignedInteger> NttRlev<T> {
    /// Creates a new [`NttRlev<T>`].
    #[inline]
    pub fn new(data: Vec<NttRlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`NttRlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, poly_length: usize) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| NttRlwe::zero(poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|rlwe| rlwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`NttRlev<T>`].
    #[inline]
    pub fn data(&self) -> &[NttRlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`NttRlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, NttRlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`NttRlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, NttRlwe<T>> {
        self.data.iter_mut()
    }
}
