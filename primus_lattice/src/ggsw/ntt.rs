use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::Ggsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
#[derive(Clone)]
pub struct NttGgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(NttGgsw<S, T>);
impl_bytes_conversion!(NttGgsw<S, T>);
impl_zero!(NttGgsw<S, T>);
impl_basic_operation_single_modulus!(NttGgsw<S, T>);
impl_intt!(NttGgsw<S, T>, Ggsw);
