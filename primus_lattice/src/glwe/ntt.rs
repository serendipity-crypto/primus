use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use super::Glwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
#[derive(Clone)]
pub struct NttGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(NttGlwe<S, T>);
impl_bytes_conversion!(NttGlwe<S, T>);
impl_zero!(NttGlwe<S, T>);
impl_basic_operation_single_modulus!(NttGlwe<S, T>);
impl_intt!(NttGlwe<S, T>, Glwe);

impl<S, T> NttGlwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> NttGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs a modular multiplication on the `self` [`NttGlwe<S>`] with another `polynomial` [`NttPolynomial<A>`].
    #[inline]
    pub fn mul_ntt_polynomial_assign<M, A>(&mut self, polynomial: &NttPolynomial<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let poly_length = polynomial.poly_length();

        self.data.chunks_exact_mut(poly_length).for_each(|p| {
            NttPolynomial(ArrayBase(p)).mul_assign(polynomial, modulus);
        });
    }
}

impl<S, T> NttGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`NttGlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Performs a modular multiplication on the `self` [`NttGlwe<S>`] with another `polynomial` [`NttPolynomial`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, A, B>(
        &self,
        polynomial: &NttPolynomial<A>,
        result: &mut NttGlwe<B>,
        modulus: M,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = polynomial.poly_length();

        self.data
            .chunks_exact(poly_length)
            .zip(result.data.chunks_exact_mut(poly_length))
            .for_each(|(x, y)| {
                NttPolynomial(ArrayBase(x)).mul_inplace(
                    polynomial,
                    &mut NttPolynomial(ArrayBase(y)),
                    modulus,
                );
            });
    }
}
