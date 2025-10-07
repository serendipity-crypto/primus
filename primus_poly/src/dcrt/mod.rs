use num_traits::ConstZero;
use primus_integer::UnsignedInteger;
use primus_utils::Size;
use serde::{Deserialize, Serialize};

use crate::NttPolynomial;

mod basic;

mod add;
mod inv;
mod mul;
mod neg;
mod sub;

/// This structure is used to store polynomials
/// with large integer coefficients and speed up multiplication.
///
/// By the Chinese remainder theorem, a large integer
/// can be decomposed into several remainders.
///
/// If all the coefficients of a polynomial are decomposed in the same way,
/// several polynomials with relatively small coefficients can be obtained,
/// and the latter has better performance in addition and subtraction computation.
///
/// Also, applying number theory transform to each factorized polynomial,
/// we can get polynomials that are more efficient in addition, subtraction and multiplication.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DcrtPolynomial<T> {
    polys: Vec<NttPolynomial<T>>,
}

impl<T> DcrtPolynomial<T> {
    /// Creates a new [`DcrtPolynomial<T>`].
    #[inline]
    pub fn new(polys: Vec<NttPolynomial<T>>) -> Self {
        Self { polys }
    }

    /// Returns an iterator that allows reading each rns polynomial of the double crt polynomial.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, NttPolynomial<T>> {
        self.polys.iter()
    }

    /// Returns an iterator that allows modifying each rns polynomial of the double crt polynomial.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, NttPolynomial<T>> {
        self.polys.iter_mut()
    }

    /// .
    #[inline]
    pub const fn as_slice(&self) -> &[NttPolynomial<T>] {
        self.polys.as_slice()
    }

    /// .
    #[inline]
    pub const fn as_mut_slice(&mut self) -> &mut [NttPolynomial<T>] {
        self.polys.as_mut_slice()
    }
}

impl<T> DcrtPolynomial<T>
where
    T: Copy + ConstZero,
{
    /// Creates a [`DcrtPolynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(moduli_count: usize, poly_length: usize) -> Self {
        Self {
            polys: (0..moduli_count)
                .map(|_| NttPolynomial::zero(poly_length))
                .collect(),
        }
    }

    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.polys.iter().all(NttPolynomial::is_zero)
    }

    /// Sets `self` to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.polys.iter_mut().for_each(NttPolynomial::set_zero);
    }
}

impl<T: Copy> DcrtPolynomial<T> {
    /// Copy the coefficients from another slice.
    #[inline]
    pub fn copy_from<I: AsRef<[T]>>(&mut self, src: impl IntoIterator<Item = I>) {
        self.polys
            .iter_mut()
            .zip(src)
            .for_each(|(a, b)| a.copy_from(b.as_ref()));
    }
}

impl<T: UnsignedInteger> Size for DcrtPolynomial<T> {
    #[inline]
    fn size(&self) -> usize {
        self.polys.len() * self.polys[0].size()
    }
}
