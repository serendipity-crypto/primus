use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

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
