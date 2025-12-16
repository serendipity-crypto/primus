use primus_integer::{Data, DataMut, DataOwned, RawData, UnsignedInteger};
use primus_ntt::NttTable;
use primus_poly::ArrayBase;
use primus_reduce::FieldContext;

use crate::glwe::{GlweIter, GlweIterMut};

use super::NttGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::glwe::Glwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct Glev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(Glev<S, T>);
impl_bytes_conversion!(Glev<S, T>);
impl_zero!(Glev<S, T>);
impl_iters!(Glev);
impl_iter_sub_structure!(Glev<S, T>, Glwe);
impl_basic_operation_single_modulus!(Glev<S, T>);
impl_ntt!(Glev<S, T>, NttGlev);
