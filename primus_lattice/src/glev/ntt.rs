use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::Glev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct NttGlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(NttGlev<S, T>);
impl_bytes_conversion!(NttGlev<S, T>);
impl_zero!(NttGlev<S, T>);
impl_basic_operation_single_modulus!(NttGlev<S, T>);
impl_intt!(NttGlev<S, T>, Glev);
