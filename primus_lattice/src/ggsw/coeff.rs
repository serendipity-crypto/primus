use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::NttGgsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::glev::Glev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct Ggsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Ggsw<S, T>);
impl_bytes_conversion!(Ggsw<S, T>);
impl_zero!(Ggsw<S, T>);
impl_iter_sub_structure!(Ggsw<S, T>, glev);
impl_basic_operation_single_modulus!(Ggsw<S, T>);
impl_ntt!(Ggsw<S, T>, NttGgsw);
