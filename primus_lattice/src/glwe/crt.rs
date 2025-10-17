use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use super::DcrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
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
impl_basic_operation_multiple_modulus!(CrtGlwe<S, T>);
impl_crt_ntt!(CrtGlwe<S, T>, DcrtGlwe);

impl<S, T> CrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> CrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> CrtGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`CrtGlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtGlwe<T>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        moduli: &[M],
        table: &Table,
        cipher_single_modulus_len: usize,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        izip!(
            result.data.chunks_exact_mut(cipher_single_modulus_len),
            dcrt_polynomial.iter_each_modulus(poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, poly, ntt_table, modulus)| {
            glwe.chunks_exact_mut(poly_length).for_each(|a| {
                ntt_table.transform_slice(a);
                NttPolynomial(ArrayBase(a)).mul_assign(&NttPolynomial(ArrayBase(poly)), *modulus);
            });
        });
    }
}
