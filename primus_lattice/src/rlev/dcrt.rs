use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::DcrtRlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct DcrtRlev<T: UnsignedInteger> {
    data: Vec<DcrtRlwe<T>>,
}

impl<T: UnsignedInteger> DcrtRlev<T> {
    /// Creates a new [`DcrtRlev<T>`].
    #[inline]
    pub fn new(data: Vec<DcrtRlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`DcrtRlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, moduli_count: usize, poly_length: usize) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| DcrtRlwe::zero(moduli_count, poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|rlwe| rlwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`DcrtRlev<T>`].
    #[inline]
    pub fn data(&self) -> &[DcrtRlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`DcrtRlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, DcrtRlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`DcrtRlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, DcrtRlwe<T>> {
        self.data.iter_mut()
    }
}
