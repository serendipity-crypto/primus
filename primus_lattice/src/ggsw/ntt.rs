use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::glev::{NttGlevIter, NttGlevIterMut};

use super::Ggsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::glev::NttGlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct NttGgsw<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(NttGgsw<S, T>);
impl_bytes_conversion!(NttGgsw<S, T>);
impl_zero!(NttGgsw<S, T>);
impl_iters!(NttGgsw);
impl_iter_sub_structure!(NttGgsw<S, T>, NttGlev);
impl_basic_operation_single_modulus!(NttGgsw<S, T>);
impl_intt!(NttGgsw<S, T>, Ggsw);
