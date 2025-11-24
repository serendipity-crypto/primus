use primus_integer::UnsignedInteger;
use primus_ntt::NttTable;
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlwe::{RlweIter, RlweIterMut};

use super::NttRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::rlwe::Rlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct Rlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(Rlev<S, T>);
impl_bytes_conversion!(Rlev<S, T>);
impl_zero!(Rlev<S, T>);
impl_iters!(Rlev);
impl_iter_sub_structure!(Rlev<S, T>, Rlwe);
impl_basic_operation_single_modulus!(Rlev<S, T>);
impl_ntt!(Rlev<S, T>, NttRlev);
