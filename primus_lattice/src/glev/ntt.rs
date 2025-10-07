use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::NttGlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different powers
/// of a base, used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct NttGlev<T: UnsignedInteger> {
    /// A vector of RLWE ciphertexts, each encrypted message with a different `basis`.
    data: Vec<NttGlwe<T>>,
}

impl<T: UnsignedInteger> NttGlev<T> {
    /// Creates a new [`NttGlev<T>`].
    #[inline]
    pub fn new(data: Vec<NttGlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`NttGlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, dimension: usize, poly_length: usize) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| NttGlwe::zero(dimension, poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|glwe| glwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`NttGlev<T>`].
    #[inline]
    pub fn data(&self) -> &[NttGlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`NttGlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, NttGlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`NttGlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, NttGlwe<T>> {
        self.data.iter_mut()
    }
}
