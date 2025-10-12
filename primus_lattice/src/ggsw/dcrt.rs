use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::CrtGgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone)]
pub struct DcrtGgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtGgsw<S, T>);
impl_bytes_conversion!(DcrtGgsw<S, T>);
impl_zero!(DcrtGgsw<S, T>);
impl_basic_operation_multiple_modulus!(DcrtGgsw<S, T>);
impl_crt_intt!(DcrtGgsw<S, T>, CrtGgsw);
