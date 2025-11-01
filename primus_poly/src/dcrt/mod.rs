use primus_integer::{ByteCount, UnsignedInteger, izip, size::Size};
use primus_reduce::ops::ReduceMulAdd;

use crate::{ArrayBase, Data, DataMut, DataOwned, RawData};

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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DcrtPolynomial<S, T = <S as RawData>::Elem>(pub ArrayBase<S, T>)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`DcrtPolynomial<T>`].
    #[inline]
    pub fn new(polys: ArrayBase<S, T>) -> Self {
        Self(polys)
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`DcrtPolynomial<T>`] with all coefficients equal to zero.
    #[inline]
    pub fn zero(dcrt_poly_length: usize) -> Self {
        Self(ArrayBase::zero(dcrt_poly_length))
    }

    #[inline]
    pub fn into_owned(self) -> S {
        self.0.0
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Returns an iterator that allows modifying each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_each_modulus_mut(
        &mut self,
        poly_length: usize,
    ) -> std::slice::ChunksExactMut<'_, T> {
        self.0.chunks_exact_mut(poly_length)
    }

    /// Sets `self` to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.0.set_zero();
    }

    /// Copy the coefficients from another slice.
    #[inline]
    pub fn copy_from<A>(&mut self, src: &DcrtPolynomial<A, T>)
    where
        A: RawData<Elem = T> + Data,
    {
        self.0.copy_from_slice(src.0.as_ref());
    }

    /// Performs `self *= rhs` according to `moduli`.
    #[inline]
    pub fn add_mul_assign<M, A, B>(
        &mut self,
        b: &DcrtPolynomial<A, T>,
        c: &DcrtPolynomial<B, T>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        izip!(
            self.iter_each_modulus_mut(poly_length),
            b.iter_each_modulus(poly_length),
            c.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(xs, ys, zs, &modulus)| {
            ArrayBase(xs).add_mul_element_wise_assign(&ArrayBase(ys), &ArrayBase(zs), modulus);
        })
    }
}

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Returns an iterator that allows reading each value or coefficient of the polynomial.
    #[inline]
    pub fn iter_each_modulus(&self, poly_length: usize) -> std::slice::ChunksExact<'_, T> {
        self.0.chunks_exact(poly_length)
    }

    /// Returns `true` if `self` is equal to `0`.
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline]
    pub fn dcrt_poly_length(&self) -> usize {
        self.0.len()
    }
}

impl<S, T> Size for DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn bytes_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<S, T> AsRef<[T]> for DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<S, T> AsMut<[T]> for DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}
