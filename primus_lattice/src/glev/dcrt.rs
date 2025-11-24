use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::DcrtTable;
use primus_poly::{
    ArrayBase, BigUintPolynomial, CrtPolynomial, Data, DataMut, DataOwned, DcrtPolynomial, RawData,
};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{
    context::DcrtGlevContext,
    glwe::{DcrtGlwe, DcrtGlweIter, DcrtGlweIterMut},
};

use super::CrtGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::glwe::DcrtGlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct DcrtGlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(DcrtGlev<S, T>);
impl_bytes_conversion!(DcrtGlev<S, T>);
impl_zero!(DcrtGlev<S, T>);
impl_iters!(DcrtGlev);
impl_iter_sub_structure!(DcrtGlev<S, T>, DcrtGlwe);
impl_basic_operation_multiple_modulus!(DcrtGlev<S, T>);
impl_crt_intt!(DcrtGlev<S, T>, CrtGlev);

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    pub fn mul_crt_poly_inplace<M, Table, A, B>(
        &self,
        crt_poly: &CrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        basis: &BigUintApproxSignedBasis<T>,
        table: &Table,
        rns_base: &RNSBase<T, M>,
        context: &mut DcrtGlevContext<T>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = table.poly_length();
        let big_uint_value_len = rns_base.big_uint_value_len();
        let basis_value = basis.basis_value();
        let moduli = rns_base.moduli();
        let dcrt_glwe_len = result.0.len();

        let (adjust_big_uint_values, decomposed_unsigned_values, carries, multi_residues) =
            context.as_mut();

        rns_base.compose_multiple_values_inplace(
            crt_poly.as_ref(),
            adjust_big_uint_values,
            poly_length,
        );

        basis.init_value_carry_slice_inplace(adjust_big_uint_values, carries, big_uint_value_len);

        result.set_zero();

        self.iter_dcrt_glwe(dcrt_glwe_len)
            .zip(basis.decomposer_iter())
            .for_each(|(dcrt_glwe, once_decomposer)| {
                once_decomposer.unsigned_decompose_slice_inplace(
                    adjust_big_uint_values,
                    decomposed_unsigned_values,
                    carries,
                    big_uint_value_len,
                );

                rns_base.wrapping_decompose_small_values_inplace(
                    decomposed_unsigned_values,
                    multi_residues,
                    poly_length,
                    basis_value,
                );

                table.transform_slice(multi_residues);

                result.add_dcrt_glwe_mul_dcrt_polynomial_assign(
                    &dcrt_glwe,
                    &DcrtPolynomial(multi_residues.as_ref()),
                    poly_length,
                    moduli,
                );
            });
    }

    pub fn mul_big_uint_poly_inplace<M, Table, A, B>(
        &self,
        big_uint_poly: &BigUintPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        basis: &BigUintApproxSignedBasis<T>,
        table: &Table,
        rns_base: &RNSBase<T, M>,
        context: &mut DcrtGlevContext<T>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = table.poly_length();
        let big_uint_value_len = rns_base.big_uint_value_len();
        let dcrt_glwe_len = result.0.len();
        let basis_value = basis.basis_value();
        let moduli = rns_base.moduli();

        let (adjust_big_uint_values, decomposed_unsigned_values, carries, multi_residues) =
            context.as_mut();

        basis.init_value_carry_slice(
            big_uint_poly.as_slice(),
            adjust_big_uint_values,
            carries,
            big_uint_value_len,
        );

        result.set_zero();

        self.iter_dcrt_glwe(dcrt_glwe_len)
            .zip(basis.decomposer_iter())
            .for_each(|(dcrt_glwe, once_decomposer)| {
                once_decomposer.unsigned_decompose_slice_inplace(
                    adjust_big_uint_values,
                    decomposed_unsigned_values,
                    carries,
                    big_uint_value_len,
                );

                rns_base.wrapping_decompose_small_values_inplace(
                    decomposed_unsigned_values,
                    multi_residues,
                    poly_length,
                    basis_value,
                );

                table.transform_slice(multi_residues);

                result.add_dcrt_glwe_mul_dcrt_polynomial_assign(
                    &dcrt_glwe,
                    &DcrtPolynomial(multi_residues.as_ref()),
                    poly_length,
                    moduli,
                );
            });
    }
}
