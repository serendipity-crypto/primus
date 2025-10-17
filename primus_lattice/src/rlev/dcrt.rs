use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::CrtRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct DcrtRlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtRlev<S, T>);
impl_bytes_conversion!(DcrtRlev<S, T>);
impl_zero!(DcrtRlev<S, T>);
impl_basic_operation_multiple_modulus!(DcrtRlev<S, T>);
impl_crt_intt!(DcrtRlev<S, T>, CrtRlev);

impl<S, T> DcrtRlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtRlev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtRlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
}
