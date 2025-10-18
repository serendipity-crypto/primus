use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::Rgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::NttRlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct NttRgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(NttRgsw<S, T>);
impl_bytes_conversion!(NttRgsw<S, T>);
impl_zero!(NttRgsw<S, T>);
impl_iter_sub_structure!(NttRgsw<S, T>, ntt_glev);
impl_basic_operation_single_modulus!(NttRgsw<S, T>);
impl_intt!(NttRgsw<S, T>, Rgsw);
