use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::DcrtGlwe;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct DcrtGlev<T: UnsignedInteger> {
    data: Vec<DcrtGlwe<T>>,
}

impl<T: UnsignedInteger> DcrtGlev<T> {
    /// Creates a new [`DcrtGlev<T>`].
    #[inline]
    pub fn new(data: Vec<DcrtGlwe<T>>) -> Self {
        Self { data }
    }

    /// Creates a [`DcrtGlev<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(
        decompose_length: usize,
        dimension: usize,
        moduli_count: usize,
        poly_length: usize,
    ) -> Self {
        Self {
            data: (0..decompose_length)
                .map(|_| DcrtGlwe::zero(dimension, moduli_count, poly_length))
                .collect(),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.iter_mut().for_each(|glwe| glwe.set_zero());
    }

    /// Returns a reference to the `data` of this [`DcrtGlev<T>`].
    #[inline]
    pub fn data(&self) -> &[DcrtGlwe<T>] {
        self.data.as_ref()
    }

    /// Returns an iterator over the `data` of this [`DcrtGlev<T>`].
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, DcrtGlwe<T>> {
        self.data.iter()
    }

    /// Returns an iterator over the `data` of this [`DcrtGlev<T>`]
    /// that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, DcrtGlwe<T>> {
        self.data.iter_mut()
    }
}
