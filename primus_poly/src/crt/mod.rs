use num_traits::ConstZero;
use primus_integer::UnsignedInteger;
use primus_utils::Size;
use serde::{Deserialize, Serialize};

use crate::Polynomial;

mod basic;

mod add;
mod mul;
mod neg;
mod sub;

/// This structure is used to store polynomials with large integer coefficients.
///
/// By the Chinese remainder theorem,
/// a large integer can be decomposed into several remainders.
///
/// If all the coefficients of a polynomial are decomposed in the same way,
/// several polynomials with relatively small coefficients can be obtained,
/// and the latter has better performance in addition and subtraction computation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CrtPolynomial<T> {
    polys: Vec<Polynomial<T>>,
}

impl<T> CrtPolynomial<T> {
    /// Creates a new [`CrtPolynomial<T>`].
    #[inline]
    pub fn new(polys: Vec<Polynomial<T>>) -> Self {
        Self { polys }
    }

    #[inline]
    pub fn into_vec(self) -> Vec<Polynomial<T>> {
        self.polys
    }

    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, Polynomial<T>> {
        self.polys.iter()
    }

    /// Returns an iterator that allows modifying each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Polynomial<T>> {
        self.polys.iter_mut()
    }
}

impl<T> CrtPolynomial<T>
where
    T: Copy + ConstZero,
{
    /// Creates a [`CrtPolynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(moduli_count: usize, poly_length: usize) -> Self {
        Self {
            polys: (0..moduli_count)
                .map(|_| Polynomial::zero(poly_length))
                .collect(),
        }
    }

    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.polys.iter().all(Polynomial::is_zero)
    }

    /// Sets `self` to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.polys.iter_mut().for_each(Polynomial::set_zero);
    }
}

impl<T: Copy> CrtPolynomial<T> {
    /// Copy the coefficients from another slice.
    #[inline]
    pub fn copy_from<I: AsRef<[T]>>(&mut self, src: impl IntoIterator<Item = I>) {
        self.polys
            .iter_mut()
            .zip(src)
            .for_each(|(a, b)| a.copy_from(b.as_ref()));
    }
}

impl<T: UnsignedInteger> Size for CrtPolynomial<T> {
    #[inline]
    fn size(&self) -> usize {
        self.polys.len() * self.polys[0].size()
    }
}
