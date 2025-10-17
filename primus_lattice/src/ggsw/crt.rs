use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::DcrtGgsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
#[derive(Clone)]
pub struct CrtGgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtGgsw<S, T>);
impl_bytes_conversion!(CrtGgsw<S, T>);
impl_zero!(CrtGgsw<S, T>);
impl_basic_operation_multiple_modulus!(CrtGgsw<S, T>);
impl_crt_ntt!(CrtGgsw<S, T>, DcrtGgsw);
