use primus_integer::{Data, DataMut, DataOwned, RawData, UnsignedInteger, izip};
use primus_ntt::DcrtTable;
use primus_poly::ArrayBase;
use primus_reduce::FieldContext;

use crate::glwe::{CrtGlweIter, CrtGlweIterMut};

use super::DcrtGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--cd--|
///
/// where `c1` to `cd` are [`crate::glwe::CrtGlwe`] with same parameter, `d` is the decompose length.
#[derive(Clone)]
pub struct CrtGlev<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(CrtGlev<S, T>);
impl_bytes_conversion!(CrtGlev<S, T>);
impl_zero!(CrtGlev<S, T>);
impl_iters!(CrtGlev);
impl_iter_sub_structure!(CrtGlev<S, T>, CrtGlwe);
impl_basic_operation_multiple_modulus!(CrtGlev<S, T>);
impl_crt_ntt!(CrtGlev<S, T>, DcrtGlev);
