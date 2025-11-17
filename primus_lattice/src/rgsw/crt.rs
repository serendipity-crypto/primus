use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlev::{CrtRlevIter, CrtRlevIterMut};

use super::DcrtRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::CrtRlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct CrtRgsw<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(CrtRgsw<S, T>);
impl_bytes_conversion!(CrtRgsw<S, T>);
impl_zero!(CrtRgsw<S, T>);
impl_iters!(CrtRgsw);
impl_iter_sub_structure!(CrtRgsw<S, T>, CrtRlev);
impl_basic_operation_multiple_modulus!(CrtRgsw<S, T>);
impl_crt_ntt!(CrtRgsw<S, T>, DcrtRgsw);
