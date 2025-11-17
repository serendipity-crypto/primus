use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::glwe::{NttGlweIter, NttGlweIterMut};

use super::Glev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::glwe::NttGlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct NttGlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(NttGlev<S, T>);
impl_bytes_conversion!(NttGlev<S, T>);
impl_zero!(NttGlev<S, T>);
impl_iters!(NttGlev);
impl_iter_sub_structure!(NttGlev<S, T>, NttGlwe);
impl_basic_operation_single_modulus!(NttGlev<S, T>);
impl_intt!(NttGlev<S, T>, Glev);
