use primus_integer::UnsignedInteger;

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

mod add;
mod neg;
mod sub;

/// Represents a polynomial where coefficients are elements of a specified numeric `T`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BigUintPolynomial<S, T = <S as RawData>::Elem>(pub ArrayBase<S, T>)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`BigUintPolynomial<T>`].
    #[inline]
    pub fn new(big_uint_poly: ArrayBase<S, T>) -> Self {
        Self(big_uint_poly)
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`BigUintPolynomial<S>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(poly_length: usize, value_len: usize) -> Self {
        Self(ArrayBase::zero(poly_length * value_len))
    }

    /// Drop self, and return the vector.
    #[inline]
    pub fn into_owned(self) -> S {
        self.0.0
    }

    /// Constructs a new polynomial from a slice.
    #[inline]
    pub fn from_slice(polynomial: &[T]) -> Self {
        Self::new(ArrayBase::from_slice(polynomial))
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Copy the coefficients from another slice.
    #[inline]
    pub fn copy_from(&mut self, src: impl AsRef<[T]>) {
        self.0.copy_from_slice(src.as_ref())
    }

    /// Sets `self` to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.0.set_zero();
    }

    /// Extracts a mutable slice of the entire vector.
    ///
    /// Equivalent to `&mut s[..]`.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        self.0.as_mut()
    }

    /// Returns an iterator that allows modifying each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_mut(&mut self, value_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.0.chunks_exact_mut(value_len)
    }
}

impl<S, T> BigUintPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Extracts a slice containing the entire vector.
    ///
    /// Equivalent to `&s[..]`.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        self.0.as_ref()
    }

    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn iter(&self, value_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.0.chunks_exact(value_len)
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }
}
