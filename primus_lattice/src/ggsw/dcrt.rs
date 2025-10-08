use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

use crate::DcrtGlev;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct DcrtGgsw<T: UnsignedInteger> {
    a: Vec<DcrtGlev<T>>,
    b: DcrtGlev<T>,
}

impl<T: UnsignedInteger> DcrtGgsw<T> {
    /// Creates a new [`DcrtGgsw<T>`].
    #[inline]
    pub fn new(a: Vec<DcrtGlev<T>>, b: DcrtGlev<T>) -> Self {
        Self { a, b }
    }

    /// Creates a [`DcrtGgsw<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(
        decompose_length: usize,
        dimension: usize,
        moduli_count: usize,
        poly_length: usize,
    ) -> Self {
        Self {
            a: (0..dimension)
                .map(|_| DcrtGlev::zero(decompose_length, dimension, moduli_count, poly_length))
                .collect(),
            b: DcrtGlev::zero(decompose_length, dimension, moduli_count, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.iter_mut().for_each(|glev| glev.set_zero());
        self.b.set_zero();
    }

    pub fn a(&self) -> &[DcrtGlev<T>] {
        &self.a
    }

    pub fn a_mut(&mut self) -> &mut [DcrtGlev<T>] {
        &mut self.a
    }

    pub fn b(&self) -> &DcrtGlev<T> {
        &self.b
    }

    pub fn b_mut(&mut self) -> &mut DcrtGlev<T> {
        &mut self.b
    }
}
