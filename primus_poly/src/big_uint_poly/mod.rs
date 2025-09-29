use num_traits::{ConstZero, Zero};
use primus_integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

mod add;
mod neg;
mod sub;

/// Represents a polynomial where coefficients are elements of a specified numeric `T`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct BigUintPolynomial<T: UnsignedInteger> {
    big_uint_poly: Vec<T>,
}

impl<T: UnsignedInteger> Default for BigUintPolynomial<T> {
    #[inline]
    fn default() -> Self {
        Self {
            big_uint_poly: Vec::new(),
        }
    }
}

impl<T: UnsignedInteger> BigUintPolynomial<T> {
    /// Creates a new [`BigUintPolynomial<T>`].
    #[inline]
    pub fn new(big_uint_poly: Vec<T>) -> Self {
        Self { big_uint_poly }
    }

    /// Drop self, and return the vector.
    #[inline]
    pub fn into_vec(self) -> Vec<T> {
        self.big_uint_poly
    }

    /// Constructs a new polynomial from a slice.
    #[inline]
    pub fn from_slice(big_uint_polynomial: &[T]) -> Self {
        Self::new(big_uint_polynomial.to_vec())
    }

    /// Copy the coefficients from another slice.
    #[inline]
    pub fn copy_from(&mut self, src: impl AsRef<[T]>) {
        self.big_uint_poly.copy_from_slice(src.as_ref())
    }

    /// Creates a [`BigUintPolynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(poly_len: usize, value_len: usize) -> Self {
        Self {
            big_uint_poly: vec![<T as ConstZero>::ZERO; poly_len * value_len],
        }
    }

    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.big_uint_poly.iter().all(<T as Zero>::is_zero)
    }

    /// Sets `self` to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.big_uint_poly.fill(<T as ConstZero>::ZERO);
    }

    /// Extracts a slice containing the entire vector.
    ///
    /// Equivalent to `&s[..]`.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        self.big_uint_poly.as_slice()
    }

    /// Extracts a mutable slice of the entire vector.
    ///
    /// Equivalent to `&mut s[..]`.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        self.big_uint_poly.as_mut_slice()
    }

    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn iter(&self, value_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.big_uint_poly.chunks_exact(value_len)
    }

    /// Returns an iterator that allows modifying each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_mut(&mut self, value_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.big_uint_poly.chunks_exact_mut(value_len)
    }
}
