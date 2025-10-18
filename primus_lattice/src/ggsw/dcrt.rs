use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::CrtGgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::glev::DcrtGlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct DcrtGgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtGgsw<S, T>);
impl_bytes_conversion!(DcrtGgsw<S, T>);
impl_zero!(DcrtGgsw<S, T>);
impl_iter_sub_structure!(DcrtGgsw<S, T>, dcrt_glev);
impl_basic_operation_multiple_modulus!(DcrtGgsw<S, T>);
impl_crt_intt!(DcrtGgsw<S, T>, CrtGgsw);
