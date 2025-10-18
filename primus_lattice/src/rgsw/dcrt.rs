use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::CrtRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// ## Structure of the `data`
///
/// |--c1--|....|--ck--|--c[k+1]--|
///
/// where `c1` to `c[k+1]` are [`crate::rlev::DcrtRlev`] with same parameter, `k` is the dimension.
#[derive(Clone)]
pub struct DcrtRgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtRgsw<S, T>);
impl_bytes_conversion!(DcrtRgsw<S, T>);
impl_zero!(DcrtRgsw<S, T>);
impl_iter_sub_structure!(DcrtRgsw<S, T>, dcrt_glev);
impl_basic_operation_multiple_modulus!(DcrtRgsw<S, T>);
impl_crt_intt!(DcrtRgsw<S, T>, CrtRgsw);

impl<S, T> DcrtRgsw<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_dcrt_rgsw_mut(
        &mut self,
        dcrt_rgsw_len: usize,
    ) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(dcrt_rgsw_len)
    }
}

impl<S, T> DcrtRgsw<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_dcrt_rgsw(&self, dcrt_rgsw_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(dcrt_rgsw_len)
    }
}
