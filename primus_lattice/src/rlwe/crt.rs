use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use super::DcrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
#[derive(Clone)]
pub struct CrtRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtRlwe<S, T>);
impl_bytes_conversion!(CrtRlwe<S, T>);
impl_zero!(CrtRlwe<S, T>);
impl_basic_operation_multiple_modulus!(CrtRlwe<S, T>);
impl_crt_ntt!(CrtRlwe<S, T>, DcrtRlwe);

impl<S, T> CrtRlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`CrtRlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtRlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtRlwe<B>,
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
            .data
            .chunks_exact_mut(crt_poly_length)
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
