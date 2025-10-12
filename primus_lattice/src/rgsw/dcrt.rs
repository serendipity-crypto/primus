use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::CrtRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone)]
pub struct DcrtRgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtRgsw<S, T>);
impl_bytes_conversion!(DcrtRgsw<S, T>);
impl_zero!(DcrtRgsw<S, T>);
impl_basic_operation_multiple_modulus!(DcrtRgsw<S, T>);
impl_crt_intt!(DcrtRgsw<S, T>, CrtRgsw);
