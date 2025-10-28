use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_factor::ShoupFactor;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{context::DcrtGlevContext, ggsw::DcrtGgsw, glev::DcrtGlev};

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

    pub fn mul_scalar_assign<M>(
        &mut self,
        scalar_residue: &[T],
        poly_length: usize,
        crt_poly_len: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
    {
        self.iter_crt_poly_mut(crt_poly_len).for_each(|crt_poly| {
            CrtPolynomial(ArrayBase(crt_poly)).mul_scalar_assign(
                scalar_residue,
                poly_length,
                moduli,
            );
        });
    }

    /// Perform `self = self * X^r`.
    pub fn mul_monic_monomial_assign<M>(
        &mut self,
        r: usize,
        poly_length: usize,
        crt_poly_len: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
    {
        if r < poly_length {
            let n_sub_r = poly_length - r;
            let rotate = |poly: &mut [T], modulus: M| {
                poly.rotate_right(r);
                poly[0..n_sub_r]
                    .iter_mut()
                    .for_each(|v| modulus.reduce_neg_assign(v));
            };

            self.iter_crt_poly_mut(crt_poly_len).for_each(|crt_poly| {
                crt_poly
                    .chunks_exact_mut(poly_length)
                    .zip(moduli)
                    .for_each(|(poly, &modulus)| rotate(poly, modulus));
            });
        } else {
            let r = (poly_length << 1) - r;
            let n_sub_r = poly_length.checked_sub(r).expect("r > 2N !");

            let rotate = |poly: &mut [T], modulus: M| {
                poly.rotate_left(r);
                poly[n_sub_r..]
                    .iter_mut()
                    .for_each(|v| modulus.reduce_neg_assign(v));
            };

            self.iter_crt_poly_mut(crt_poly_len).for_each(|crt_poly| {
                crt_poly
                    .chunks_exact_mut(poly_length)
                    .zip(moduli)
                    .for_each(|(poly, &modulus)| rotate(poly, modulus));
            });
        }
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

    pub fn mul_scalar_inplace<M, A>(
        &self,
        scalar_residue: &[T],
        result: &mut CrtGlwe<A>,
        poly_length: usize,
        crt_poly_len: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.iter_crt_poly(crt_poly_len)
            .zip(result.iter_crt_poly_mut(crt_poly_len))
            .for_each(|(in_crt_poly, out_crt_poly)| {
                CrtPolynomial(ArrayBase(in_crt_poly)).mul_scalar_inplace(
                    scalar_residue,
                    &mut CrtPolynomial(ArrayBase(out_crt_poly)),
                    poly_length,
                    moduli,
                );
            });
    }

    pub fn mul_factor_inplace<A>(
        &self,
        scalar: &[ShoupFactor<T>],
        result: &mut CrtGlwe<A>,
        poly_length: usize,
        crt_poly_len: usize,
        moduli: &[T],
    ) where
        A: RawData<Elem = T> + DataMut,
    {
        self.iter_crt_poly(crt_poly_len)
            .zip(result.iter_crt_poly_mut(crt_poly_len))
            .for_each(|(in_crt_poly, out_crt_poly)| {
                CrtPolynomial(ArrayBase(in_crt_poly)).mul_factor_inplace(
                    scalar,
                    &mut CrtPolynomial(ArrayBase(out_crt_poly)),
                    poly_length,
                    moduli,
                );
            });
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

    pub fn mul_dcrt_ggsw_inplace<M, Table, A, B>(
        &self,
        dcrt_ggsw: &DcrtGgsw<A>,
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
        let crt_poly_len = table.crt_poly_length();
        let dcrt_glev_len = basis.decompose_length() * self.as_ref().len();

        result.set_zero();

        dcrt_ggsw
            .iter_dcrt_glev(dcrt_glev_len)
            .zip(self.iter_crt_poly(crt_poly_len))
            .for_each(|(dcrt_glev, crt_poly)| {
                result.add_dcrt_glev_mul_crt_poly_assign(
                    &DcrtGlev::new(ArrayBase(dcrt_glev)),
                    &CrtPolynomial(ArrayBase(crt_poly)),
                    basis,
                    table,
                    rns_base,
                    context,
                );
            });
    }
}
