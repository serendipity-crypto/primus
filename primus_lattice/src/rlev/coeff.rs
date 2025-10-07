use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::Rlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different powers
/// of a base, used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct Rlev<T: UnsignedInteger> {
    /// A vector of RLWE ciphertexts, each encrypted message with a different `basis`.
    data: Vec<Rlwe<T>>,
}

impl<T: UnsignedInteger> Rlev<T> {
    /// Creates a new [`Rlev<T>`].
    #[inline]
    pub fn new(data: Vec<Rlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`Rlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, poly_length: usize) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| Rlwe::zero(poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|rlwe| rlwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`Rlev<T>`].
    #[inline]
    pub fn data(&self) -> &[Rlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`Rlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, Rlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`Rlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, Rlwe<T>> {
        self.data.iter_mut()
    }
}
