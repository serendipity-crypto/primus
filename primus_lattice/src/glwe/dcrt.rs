use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{context::DcrtGlevContext, glev::DcrtGlev};

use super::CrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
///
/// ## Structure of the `data`
///
/// |--a1--|....|--ak--|--b--|
///
/// where `a1`...`ak` and `b` are [`DcrtPolynomial`] with same poly length and moduli count, `k` is the dimension.
#[derive(Clone)]
pub struct DcrtGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtGlwe<S, T>);
impl_bytes_conversion!(DcrtGlwe<S, T>);
impl_zero!(DcrtGlwe<S, T>);
impl_iter_sub_structure!(DcrtGlwe<S, T>, dcrt_poly);
impl_basic_operation_multiple_modulus!(DcrtGlwe<S, T>);
impl_crt_intt!(DcrtGlwe<S, T>, CrtGlwe);

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts mutable slice of `a` and `b` of this [`DcrtGlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        unsafe { self.data.0.split_at_mut_unchecked(mid) }
    }

    pub fn neg_assign<M>(&mut self, dcrt_poly_length: usize, poly_length: usize, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.iter_dcrt_poly_mut(dcrt_poly_length)
            .for_each(|dcrt_poly| {
                DcrtPolynomial(ArrayBase(dcrt_poly)).neg_assign(poly_length, moduli);
            });
    }

    pub fn add_dcrt_glwe_mul_dcrt_polynomial_assign<M, A, B>(
        &mut self,
        dcrt_glwe: &DcrtGlwe<A>,
        dcrt_polynomial: &DcrtPolynomial<B>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        let dcrt_poly_length = dcrt_polynomial.dcrt_poly_length();
        izip!(
            self.iter_dcrt_poly_mut(dcrt_poly_length),
            dcrt_glwe.iter_dcrt_poly(dcrt_poly_length)
        )
        .for_each(|(x, y)| {
            DcrtPolynomial(ArrayBase(x)).add_mul_assign(
                &DcrtPolynomial(ArrayBase(y)),
                dcrt_polynomial,
                poly_length,
                moduli,
            );
        });
    }

    pub fn add_dcrt_glev_mul_big_uint_poly_assign<M, Table, A, B>(
        &mut self,
        dcrt_glev: &DcrtGlev<A>,
        big_uint_poly: &BigUintPolynomial<B>,
        basis: &BigUintApproxSignedBasis<T>,
        table: &Table,
        rns_base: &RNSBase<T, M>,
        context: &mut DcrtGlevContext<T>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        let poly_length = table.poly_length();
        let big_uint_value_len = rns_base.big_uint_value_len();
        let big_uint_poly_len = big_uint_poly.len();
        let basis_value = basis.basis_value();

        debug_assert_eq!(big_uint_poly_len, big_uint_value_len * poly_length);

        let moduli = rns_base.moduli();
        let dcrt_glwe_len = self.data.len();

        let (adjust_big_uint_values, decomposed_unsigned_values, carries, multi_residues) =
            context.as_mut();

        debug_assert_eq!(adjust_big_uint_values.len(), big_uint_poly_len);
        debug_assert_eq!(decomposed_unsigned_values.len(), poly_length);
        debug_assert_eq!(carries.len(), poly_length);
        debug_assert_eq!(multi_residues.len(), poly_length * moduli.len());
        debug_assert_eq!(
            dcrt_glev.as_ref().len(),
            dcrt_glwe_len * basis.decompose_length()
        );

        basis.init_value_carry_slice(
            big_uint_poly.as_slice(),
            adjust_big_uint_values.as_mut(),
            carries.as_mut(),
            big_uint_value_len,
        );

        izip!(
            dcrt_glev.iter_dcrt_glwe(dcrt_glwe_len),
            basis.decomposer_iter()
        )
        .for_each(|(dcrt_glwe, once_decomposer)| {
            once_decomposer.unsigned_decompose_slice_inplace(
                adjust_big_uint_values.as_ref(),
                decomposed_unsigned_values.as_mut(),
                carries.as_mut(),
                big_uint_value_len,
            );

            rns_base.wrapping_decompose_small_values_inplace(
                decomposed_unsigned_values.as_ref(),
                multi_residues.as_mut(),
                poly_length,
                basis_value,
            );

            table.transform_slice(multi_residues.as_mut());

            self.add_dcrt_glwe_mul_dcrt_polynomial_assign(
                &DcrtGlwe::new(ArrayBase(dcrt_glwe)),
                &DcrtPolynomial(ArrayBase(multi_residues.as_ref())),
                poly_length,
                moduli,
            );
        });
    }
}

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`DcrtGlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Performs a multiplication on the `self` [`DcrtGlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtGlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let dcrt_poly_length = dcrt_polynomial.dcrt_poly_length();

        izip!(
            self.iter_dcrt_poly(dcrt_poly_length),
            result.iter_dcrt_poly_mut(dcrt_poly_length),
        )
        .for_each(|(a, b)| {
            DcrtPolynomial(ArrayBase(a)).mul_inplace(
                dcrt_polynomial,
                &mut DcrtPolynomial(ArrayBase(b)),
                poly_length,
                moduli,
            );
        });
    }
}
