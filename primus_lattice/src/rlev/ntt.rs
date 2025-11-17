use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::rlwe::{NttRlweIter, NttRlweIterMut};

use super::Rlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::rlwe::NttRlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct NttRlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(NttRlev<S, T>);
impl_bytes_conversion!(NttRlev<S, T>);
impl_zero!(NttRlev<S, T>);
impl_iters!(NttRlev);
impl_iter_sub_structure!(NttRlev<S, T>, NttRlwe);
impl_basic_operation_single_modulus!(NttRlev<S, T>);
impl_intt!(NttRlev<S, T>, Rlev);
