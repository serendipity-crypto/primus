use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::NttGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct Glev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Glev<S, T>);
impl_bytes_conversion!(Glev<S, T>);
impl_zero!(Glev<S, T>);
impl_basic_operation_single_modulus!(Glev<S, T>);
impl_ntt!(Glev<S, T>, NttGlev);

impl<S, T> Glev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> Glev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> Glev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
}
