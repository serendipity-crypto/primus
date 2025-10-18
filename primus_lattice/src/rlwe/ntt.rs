use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use super::Rlwe;

pub type NttRlweOwned<T> = NttRlwe<Vec<T>>;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
///
/// ## Structure of the `data`
///
/// |------a------|------b------|
///
/// where `a` and `b` are [`NttPolynomial`] with same poly length.
#[derive(Clone)]
pub struct NttRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(NttRlwe<S, T>);
impl_bytes_conversion!(NttRlwe<S, T>);
impl_zero!(NttRlwe<S, T>);
impl_iter_sub_structure!(NttRlwe<S, T>, ntt_poly);
impl_basic_operation_single_modulus!(NttRlwe<S, T>);
impl_intt!(NttRlwe<S, T>, Rlwe);

impl<S, T> NttRlwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwe<S>`] with reference of [`NttPolynomial<A>`].
    #[inline]
    pub fn from_ref<A>(a: &NttPolynomial<A>, b: &NttPolynomial<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: ArrayBase::from_vec([a.0.as_ref(), b.0.as_ref()].concat()),
        }
    }
}

impl<S, T> NttRlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts mutable slice of `a` and `b` of this [`NttRlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_mut_unchecked(mid) }
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<S>`] with another `polynomial` [`NttPolynomial<A>`].
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

impl<S, T> NttRlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`NttRlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self) -> (&[T], &[T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<S>`] with another `polynomial` [`NttPolynomial`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, A, B>(
        &self,
        polynomial: &NttPolynomial<A>,
        result: &mut NttRlwe<B>,
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
