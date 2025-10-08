use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::CrtGlev;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct CrtGgsw<T: UnsignedInteger> {
    a: Vec<CrtGlev<T>>,
    b: CrtGlev<T>,
}

impl<T: UnsignedInteger> CrtGgsw<T> {
    /// Creates a new [`CrtGgsw<T>`].
    #[inline]
    pub fn new(a: Vec<CrtGlev<T>>, b: CrtGlev<T>) -> Self {
        Self { a, b }
    }

    /// Creates a [`CrtGgsw<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(
        decompose_length: usize,
        dimension: usize,
        moduli_count: usize,
        poly_length: usize,
    ) -> Self {
        Self {
            a: (0..dimension)
                .map(|_| CrtGlev::zero(decompose_length, dimension, moduli_count, poly_length))
                .collect(),
            b: CrtGlev::zero(decompose_length, dimension, moduli_count, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.iter_mut().for_each(|glev| glev.set_zero());
        self.b.set_zero();
    }

    pub fn a(&self) -> &[CrtGlev<T>] {
        &self.a
    }

    pub fn a_mut(&mut self) -> &mut [CrtGlev<T>] {
        &mut self.a
    }

    pub fn b(&self) -> &CrtGlev<T> {
        &self.b
    }

    pub fn b_mut(&mut self) -> &mut CrtGlev<T> {
        &mut self.b
    }
}
