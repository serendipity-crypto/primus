use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlwe::{DcrtRlweIter, DcrtRlweIterMut};

use super::CrtRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::rlwe::DcrtRlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct DcrtRlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(DcrtRlev<S, T>);
impl_bytes_conversion!(DcrtRlev<S, T>);
impl_zero!(DcrtRlev<S, T>);
impl_iters!(DcrtRlev);
impl_iter_sub_structure!(DcrtRlev<S, T>, DcrtRlwe);
impl_basic_operation_multiple_modulus!(DcrtRlev<S, T>);
impl_crt_intt!(DcrtRlev<S, T>, CrtRlev);
