use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlev::{DcrtRlevIter, DcrtRlevIterMut};

use super::CrtRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::DcrtRlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct DcrtRgsw<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(DcrtRgsw<S, T>);
impl_bytes_conversion!(DcrtRgsw<S, T>);
impl_zero!(DcrtRgsw<S, T>);
impl_iters!(DcrtRgsw);
impl_iter_sub_structure!(DcrtRgsw<S, T>, DcrtRlev);
impl_basic_operation_multiple_modulus!(DcrtRgsw<S, T>);
impl_crt_intt!(DcrtRgsw<S, T>, CrtRgsw);
