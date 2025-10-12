use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, Polynomial, RawData};
use primus_reduce::FieldContext;

use crate::NttRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
#[derive(Clone)]
pub struct Rlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Rlwe<S, T>);
impl_bytes_conversion!(Rlwe<S, T>);
impl_zero!(Rlwe<S, T>);
impl_basic_operation_single_modulus!(Rlwe<S, T>);
impl_ntt!(Rlwe<S, T>, NttRlwe);

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Rlwe<S>`] with reference of [`Polynomial<A>`].
    #[inline]
    pub fn from_ref<A>(a: &Polynomial<A>, b: &Polynomial<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: ArrayBase::from_vec([a.0.as_ref(), b.0.as_ref()].concat()),
        }
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts mutable slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.0.split_at_mut_unchecked(mid) }
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self) -> (&[T], &[T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Performs a multiplication on the `self` [`Rlwe<S>`] with another `ntt_polynomial` [`NttPolynomial<A>`],
    /// store the result into `result` [`NttRlwe<B>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table, A, B>(
        &self,
        ntt_polynomial: &NttPolynomial<A>,
        result: &mut NttRlwe<B>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = ntt_table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
            NttPolynomial(ArrayBase(poly)).mul_assign(ntt_polynomial, modulus);
        });
    }
}
