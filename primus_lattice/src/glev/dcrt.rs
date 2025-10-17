use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, BigUintPolynomial, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::glwe::DcrtGlwe;

use super::CrtGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct DcrtGlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtGlev<S, T>);
impl_bytes_conversion!(DcrtGlev<S, T>);
impl_zero!(DcrtGlev<S, T>);
impl_basic_operation_multiple_modulus!(DcrtGlev<S, T>);
impl_crt_intt!(DcrtGlev<S, T>, CrtGlev);

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    pub fn mul_polynomial_inplace<M, Table, A, B>(
        &self,
        poly: &BigUintPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        basis: &BigUintApproxSignedBasis<T>,
        table: &Table,
        rns_base: &RNSBase<T, M>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        result.set_zero();

        let poly_length = table.poly_length();
        let moduli = rns_base.moduli();
        let moduli_count = moduli.len();

        let value_len = rns_base.single_value_len();

        let mut adjust_values = vec![T::ZERO; poly.0.len()];
        let mut decomposed_values = vec![T::ZERO; poly.0.len()];
        let mut carries = vec![false; poly.0.len()];
        let mut multi_residues = vec![T::ZERO; poly_length * moduli_count];

        basis.init_value_carry_slice(poly.as_slice(), &mut adjust_values, &mut carries, value_len);

        basis.decompose_iter().for_each(|once_decompose| {
            once_decompose.decompose_slice_inplace(
                &adjust_values,
                &mut decomposed_values,
                &mut carries,
                value_len,
            );

            rns_base.decompose_multiple_values_inplace(
                &decomposed_values,
                &mut multi_residues,
                poly_length,
            );
        });
    }
}

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
}
