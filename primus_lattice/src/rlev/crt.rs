use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlwe::{CrtRlweIter, CrtRlweIterMut};

use super::DcrtRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::rlwe::CrtRlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct CrtRlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(CrtRlev<S, T>);
impl_bytes_conversion!(CrtRlev<S, T>);
impl_zero!(CrtRlev<S, T>);
impl_iters!(CrtRlev);
impl_iter_sub_structure!(CrtRlev<S, T>, CrtRlwe);
impl_basic_operation_multiple_modulus!(CrtRlev<S, T>);
impl_crt_ntt!(CrtRlev<S, T>, DcrtRlev);
