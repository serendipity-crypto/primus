use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{context::DcrtGlevContext, glwe::DcrtGlwe};

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
impl_iter_sub_structure!(DcrtGlev<S, T>, dcrt_glwe);
impl_basic_operation_multiple_modulus!(DcrtGlev<S, T>);
impl_crt_intt!(DcrtGlev<S, T>, CrtGlev);

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
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
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        result.set_zero();

        let poly_length = table.poly_length();
        let big_uint_value_len = rns_base.big_uint_value_len();
        let big_uint_poly_len = big_uint_poly.len();

        debug_assert_eq!(big_uint_poly_len, big_uint_value_len * poly_length);

        let moduli = rns_base.moduli();
        let dcrt_glwe_len = result.data.len();

        let (adjust_big_uint_values, decomposed_big_uint_values, carries, multi_residues) =
            context.as_mut();

        basis.init_value_carry_slice(
            big_uint_poly.as_slice(),
            adjust_big_uint_values.as_mut(),
            carries.as_mut(),
            big_uint_value_len,
        );

        izip!(self.iter_dcrt_glwe(dcrt_glwe_len), basis.decomposer_iter()).for_each(
            |(dcrt_glwe, once_decomposer)| {
                once_decomposer.decompose_slice_inplace(
                    adjust_big_uint_values.as_ref(),
                    decomposed_big_uint_values.as_mut(),
                    carries.as_mut(),
                    big_uint_value_len,
                );

                rns_base.decompose_big_uint_values_inplace(
                    decomposed_big_uint_values.as_ref(),
                    multi_residues.as_mut(),
                    poly_length,
                );

                table.transform_slice(multi_residues.as_mut());

                result.add_dcrt_glwe_mul_dcrt_polynomial_assign(
                    &DcrtGlwe::new(ArrayBase(dcrt_glwe)),
                    &DcrtPolynomial(ArrayBase(multi_residues.as_ref())),
                    poly_length,
                    moduli,
                );
            },
        );
    }
}
