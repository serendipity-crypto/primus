use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::NttRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct Rlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Rlev<S, T>);
impl_bytes_conversion!(Rlev<S, T>);
impl_zero!(Rlev<S, T>);
impl_basic_operation_single_modulus!(Rlev<S, T>);
impl_ntt!(Rlev<S, T>, NttRlev);
