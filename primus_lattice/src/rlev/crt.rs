use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::CrtRlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct CrtRlev<T: UnsignedInteger> {
    data: Vec<CrtRlwe<T>>,
}

impl<T: UnsignedInteger> CrtRlev<T> {
    /// Creates a new [`CrtRlev<T>`].
    #[inline]
    pub fn new(data: Vec<CrtRlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`CrtRlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, moduli_count: usize, poly_length: usize) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| CrtRlwe::zero(moduli_count, poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|rlwe| rlwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`CrtRlev<T>`].
    #[inline]
    pub fn data(&self) -> &[CrtRlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`CrtRlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, CrtRlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`CrtRlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, CrtRlwe<T>> {
        self.data.iter_mut()
    }
}
