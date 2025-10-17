use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::CrtGgsw;

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
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
impl_basic_operation_multiple_modulus!(DcrtGgsw<S, T>);
impl_crt_intt!(DcrtGgsw<S, T>, CrtGgsw);

impl<S, T> DcrtGgsw<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_dcrt_glev_mut(
        &mut self,
        dcrt_glev_len: usize,
    ) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(dcrt_glev_len)
    }
}

impl<S, T> DcrtGgsw<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_dcrt_glev(&self, dcrt_glev_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(dcrt_glev_len)
    }
}
