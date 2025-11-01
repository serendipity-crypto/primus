use primus_integer::{ByteCount, UnsignedInteger, izip, size::Size};
use primus_reduce::{lazy_ops::LazyReduceMulAdd, ops::ReduceMulAdd};

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

mod basic;
mod random;

mod add;
mod inv;
mod mul;
mod neg;
mod sub;

pub type NttPolynomialOwned<T> = NttPolynomial<Vec<T>, T>;
pub type NttPolynomialRef<'a, T> = NttPolynomial<&'a [T], T>;
pub type NttPolynomialMut<'a, T> = NttPolynomial<&'a mut [T], T>;

/// Represents a ntt polynomial where values are elements of a specified numeric `T`.
/// It stores the values of the polynomial at some particular points.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NttPolynomial<S, T = <S as RawData>::Elem>(pub ArrayBase<S, T>)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`NttPolynomial<T>`].
    #[inline]
    pub fn new(values: ArrayBase<S, T>) -> Self {
        Self(values)
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`NttPolynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self(ArrayBase::zero(poly_length))
    }

    /// Drop self, and return the data.
    #[inline]
    pub fn into_owned(self) -> S {
        self.0.0
    }

    /// Constructs a new ntt polynomial from a slice.
    #[inline]
    pub fn from_slice(polynomial: &[T]) -> Self {
        Self::new(ArrayBase::from_slice(polynomial))
    }
}

impl<S, T> NttPolynomial<S, T>
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
    pub fn iter_mut(&mut self) -> core::slice::IterMut<'_, T> {
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

    /// Performs `self = self + (a * b)`.
    #[inline]
    pub fn add_mul_assign<M, A, B>(
        &mut self,
        a: &NttPolynomial<A, T>,
        b: &NttPolynomial<B, T>,
        modulus: M,
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        izip!(self.0.iter_mut(), a.0.iter(), b.0.iter())
            .for_each(|(z, &x, &y)| *z = modulus.reduce_mul_add(x, y, *z));
    }

    /// Performs `self = self + (a * b)`.
    #[inline]
    pub fn add_mul_assign_fast<M, A>(
        &mut self,
        a: &NttPolynomial<A, T>,
        b: &NttPolynomial<A, T>,
        modulus: M,
    ) where
        M: Copy + LazyReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(self.0.iter_mut(), a.0.iter(), b.0.iter())
            .for_each(|(z, &x, &y)| *z = modulus.lazy_reduce_mul_add(x, y, *z));
    }
}

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts a slice containing the entire vector.
    ///
    /// Equivalent to `&s[..]`.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        self.0.as_ref()
    }

    /// Get the `coefficient counts`/`polynomial length` of polynomial.
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.0.len()
    }

    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn iter(&self) -> core::slice::Iter<'_, T> {
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

    /// Performs `result = self * b + c`.
    #[inline]
    pub fn mul_add_inplace<M, A>(
        &self,
        b: &Self,
        c: &Self,
        result: &mut NttPolynomial<A, T>,
        modulus: M,
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(result.iter_mut(), self.iter(), b.iter(), c.iter())
            .for_each(|(d, &a, &b, &c)| *d = modulus.reduce_mul_add(a, b, c));
    }

    /// Performs `result = self * b + c`.
    #[inline]
    pub fn mul_add_inplace_fast<M, A>(
        &self,
        b: &Self,
        c: &Self,
        result: &mut NttPolynomial<A, T>,
        modulus: M,
    ) where
        M: Copy + LazyReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + DataMut,
    {
        izip!(result.iter_mut(), self.iter(), b.iter(), c.iter())
            .for_each(|(d, &a, &b, &c)| *d = modulus.lazy_reduce_mul_add(a, b, c));
    }
}

impl<S, T> Size for NttPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn bytes_count(&self) -> usize {
        self.poly_length() * <T as ByteCount>::BYTES_COUNT
    }
}
