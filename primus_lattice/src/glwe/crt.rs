use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use super::DcrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
///
/// ## Structure of the `data`
///
/// |--a1--|....|--ak--|--b--|
///
/// where `a1`...`ak` and `b` are [`primus_poly::crt::CrtPolynomial`] with same poly length and moduli count, `k` is the dimension.
#[derive(Clone)]
pub struct CrtGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtGlwe<S, T>);
impl_bytes_conversion!(CrtGlwe<S, T>);
impl_zero!(CrtGlwe<S, T>);
impl_iter_sub_structure!(CrtGlwe<S, T>, crt_poly);
impl_basic_operation_multiple_modulus!(CrtGlwe<S, T>);
impl_crt_ntt!(CrtGlwe<S, T>, DcrtGlwe);

impl<S, T> CrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts mutable slice of `a` and `b` of this [`CrtGlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        self.data.0.split_at_mut(mid)
    }
}

impl<S, T> CrtGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`CrtGlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self, mid: usize) -> (&[T], &[T]) {
        self.data.split_at(mid)
    }

    /// Performs a multiplication on the `self` [`CrtGlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtGlwe<T>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        moduli: &[M],
        table: &Table,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = table.poly_length();
        let crt_poly_length = table.crt_poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result
            .iter_dcrt_poly_mut(crt_poly_length)
            .for_each(|crt_poly| {
                table.transform_slice(crt_poly);
                DcrtPolynomial(ArrayBase(crt_poly)).mul_assign(
                    dcrt_polynomial,
                    poly_length,
                    moduli,
                );
            });
    }
}
