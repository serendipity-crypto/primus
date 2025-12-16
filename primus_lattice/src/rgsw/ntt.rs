use primus_integer::{Data, DataMut, DataOwned, RawData, UnsignedInteger};
use primus_ntt::NttTable;
use primus_poly::ArrayBase;
use primus_reduce::FieldContext;

use crate::rlev::{NttRlevIter, NttRlevIterMut};

use super::Rgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::NttRlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct NttRgsw<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(NttRgsw<S, T>);
impl_bytes_conversion!(NttRgsw<S, T>);
impl_zero!(NttRgsw<S, T>);
impl_iters!(NttRgsw);
impl_iter_sub_structure!(NttRgsw<S, T>, NttRlev);
impl_basic_operation_single_modulus!(NttRgsw<S, T>);
impl_intt!(NttRgsw<S, T>, Rgsw);
