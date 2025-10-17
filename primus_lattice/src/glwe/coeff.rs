use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use super::NttGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
#[derive(Clone)]
pub struct Glwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Glwe<S, T>);
impl_bytes_conversion!(Glwe<S, T>);
impl_zero!(Glwe<S, T>);
impl_basic_operation_single_modulus!(Glwe<S, T>);
impl_ntt!(Glwe<S, T>, NttGlwe);

impl<S, T> Glwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> Glwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> Glwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`Glwe<S>`] with another `ntt_polynomial` [`NttPolynomial<A>`],
    /// store the result into `result` [`NttGlwe<B>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table, A, B>(
        &self,
        ntt_polynomial: &NttPolynomial<A>,
        result: &mut NttGlwe<B>,
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
