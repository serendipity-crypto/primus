use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::CrtRlev;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct CrtRgsw<T: UnsignedInteger> {
    a: CrtRlev<T>,
    b: CrtRlev<T>,
}

impl<T: UnsignedInteger> CrtRgsw<T> {
    /// Creates a new [`CrtRgsw<T>`].
    #[inline]
    pub fn new(a: CrtRlev<T>, b: CrtRlev<T>) -> Self {
        Self { a, b }
    }

    /// Creates a [`CrtRgsw<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, moduli_count: usize, poly_length: usize) -> Self {
        Self {
            a: CrtRlev::zero(decompose_length, moduli_count, poly_length),
            b: CrtRlev::zero(decompose_length, moduli_count, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.set_zero();
        self.b.set_zero();
    }

    /// Returns a reference to the `-s*m` of this [`CrtRgsw<T>`].
    #[inline]
    pub fn a(&self) -> &CrtRlev<T> {
        &self.a
    }

    /// Returns a mutable reference to the `-s*m` of this [`CrtRgsw<T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut CrtRlev<T> {
        &mut self.a
    }

    /// Returns a reference to the `m` of this [`CrtRgsw<T>`].
    #[inline]
    pub fn b(&self) -> &CrtRlev<T> {
        &self.b
    }

    /// Returns a mutable reference to the `m` of this [`CrtRgsw<T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut CrtRlev<T> {
        &mut self.b
    }
}
