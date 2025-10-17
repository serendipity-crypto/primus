use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::DcrtRgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone)]
pub struct CrtRgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtRgsw<S, T>);
impl_bytes_conversion!(CrtRgsw<S, T>);
impl_zero!(CrtRgsw<S, T>);
impl_basic_operation_multiple_modulus!(CrtRgsw<S, T>);
impl_crt_ntt!(CrtRgsw<S, T>, DcrtRgsw);

impl<S, T> CrtRgsw<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_rgsw_mut(&mut self, crt_rgsw_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(crt_rgsw_len)
    }
}

impl<S, T> CrtRgsw<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_rgsw(&self, crt_rgsw_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(crt_rgsw_len)
    }
}
