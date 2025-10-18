use primus_integer::{UnsignedInteger, izip};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use super::CrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
///
/// ## Structure of the `data`
///
/// |--a1--|....|--ak--|--b--|
///
/// where `a1`...`ak` and `b` are [`primus_poly::BigUintPolynomial`] with same poly length, `k` is the dimension.
#[derive(Clone)]
pub struct BigUintGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(BigUintGlwe<S, T>);
impl_bytes_conversion!(BigUintGlwe<S, T>);
impl_zero!(BigUintGlwe<S, T>);
impl_iter_sub_structure!(BigUintGlwe<S, T>, big_uint_poly);

impl<S, T> BigUintGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn compose_assign<A, M>(
        &mut self,
        crt_glwe: &CrtGlwe<A>,
        poly_length: usize,
        crt_poly_length: usize,
        rns_base: &RNSBase<T, M>,
    ) where
        A: RawData<Elem = T> + Data,
        M: FieldContext<T>,
    {
        let single_value_len = rns_base.single_value_len();
        let big_uint_poly_length = poly_length * single_value_len;

        izip!(
            self.iter_big_uint_poly_mut(big_uint_poly_length),
            crt_glwe.iter_crt_poly(crt_poly_length),
        )
        .for_each(|(big_uint_poly, crt_poly)| {
            rns_base.compose_multiple_values_inplace(crt_poly, big_uint_poly, poly_length);
        });
    }
}

impl<S, T> BigUintGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Decomposes `self`.
    #[inline]
    pub fn decompose_inplace<A, M>(
        &self,
        result: &mut CrtGlwe<A>,
        poly_length: usize,
        crt_poly_length: usize,
        rns_base: &RNSBase<T, M>,
    ) where
        A: RawData<Elem = T> + DataMut,
        M: FieldContext<T>,
    {
        let single_value_len = rns_base.single_value_len();
        let big_uint_poly_length = poly_length * single_value_len;

        izip!(
            self.iter_big_uint_poly(big_uint_poly_length),
            result.iter_crt_poly_mut(crt_poly_length),
        )
        .for_each(|(big_uint_poly, crt_poly)| {
            rns_base.decompose_multiple_values_inplace(big_uint_poly, crt_poly, poly_length);
        });
    }
}
