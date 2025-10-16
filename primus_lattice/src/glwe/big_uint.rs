use primus_integer::UnsignedInteger;
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::CrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
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

impl<S, T> BigUintGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn compose_assign<A, M>(
        &mut self,
        crt_glwe: &CrtGlwe<A>,
        dimension: usize,
        poly_length: usize,
        rns_base: &RNSBase<T, M>,
    ) where
        A: RawData<Elem = T> + Data,
        M: FieldContext<T>,
    {
        rns_base.compose_multiple_values_inplace(
            crt_glwe.as_ref(),
            self.as_mut(),
            (dimension + 1) * poly_length,
        );
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
        dimension: usize,
        poly_length: usize,
        rns_base: &RNSBase<T, M>,
    ) where
        A: RawData<Elem = T> + DataMut,
        M: FieldContext<T>,
    {
        rns_base.decompose_multiple_values_inplace(
            self.as_ref(),
            result.as_mut(),
            (dimension + 1) * poly_length,
        );
    }
}
