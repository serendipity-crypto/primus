use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::NttRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::Rlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct Rgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Rgsw<S, T>);
impl_bytes_conversion!(Rgsw<S, T>);
impl_zero!(Rgsw<S, T>);
impl_iter_sub_structure!(Rgsw<S, T>, glev);
impl_basic_operation_single_modulus!(Rgsw<S, T>);
impl_ntt!(Rgsw<S, T>, NttRgsw);
