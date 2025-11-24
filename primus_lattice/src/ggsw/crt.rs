use primus_integer::{UnsignedInteger, izip};
use primus_ntt::DcrtTable;
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::glev::{CrtGlevIter, CrtGlevIterMut};

use super::DcrtGgsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::glev::CrtGlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct CrtGgsw<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(CrtGgsw<S, T>);
impl_bytes_conversion!(CrtGgsw<S, T>);
impl_zero!(CrtGgsw<S, T>);
impl_iters!(CrtGgsw);
impl_iter_sub_structure!(CrtGgsw<S, T>, CrtGlev);
impl_basic_operation_multiple_modulus!(CrtGgsw<S, T>);
impl_crt_ntt!(CrtGgsw<S, T>, DcrtGgsw);
