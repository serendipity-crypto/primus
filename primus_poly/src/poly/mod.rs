use num_traits::ConstZero;
use primus_integer::{UnsignedInteger, size::Size};
use primus_reduce::ops::ReduceMulAdd;

use crate::{Data, DataMut, DataOwned, RawData};

mod basic;
mod random;

mod add;
mod mul;
mod neg;
mod sub;

pub type PolynomialOwned<T> = Polynomial<Vec<T>, T>;
pub type PolynomialRef<'a, T> = Polynomial<&'a [T], T>;
pub type PolynomialMut<'a, T> = Polynomial<&'a mut [T], T>;

/// Represents a polynomial where coefficients are elements of a specified unsigned integer `T`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_iters!(Polynomial, poly);

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`Polynomial<T>`].
    #[inline]
    pub fn new(poly: S) -> Self {
        Self(poly)
    }
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`Polynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self(S::zero(poly_length))
    }

    /// Drop self, and return the vector.
    #[inline]
    pub fn into_owned(self) -> S {
        self.0
    }

    /// Constructs a new polynomial from a slice.
    #[inline]
    pub fn from_slice(polynomial: &[T]) -> Self {
        Self::new(S::from_slice(polynomial))
    }
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts a mutable slice of the entire vector.
    ///
    /// Equivalent to `&mut s[..]`.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        self.0.as_mut()
    }

    /// Returns an iterator that allows modifying each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.0.iter_mut()
    }

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
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Get the coefficient counts of polynomial.
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.0.len()
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
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn copied_iter(&self) -> core::iter::Copied<core::slice::Iter<'_, T>> {
        self.0.iter().copied()
    }

    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Evaluate the polynomial with the value `x`.
    #[inline]
    pub fn evaluate<M>(&self, x: T, modulus: M) -> T
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.0.iter().rev().fold(<T as ConstZero>::ZERO, |acc, &a| {
            modulus.reduce_mul_add(acc, x, a)
        })
    }
}

impl<S, T> Size for Polynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.byte_count()
    }
}
