use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::DcrtGgsw;

/// Represents a ciphertext in the General-GSW homomorphic encryption scheme.
#[derive(Clone)]
pub struct CrtGgsw<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtGgsw<S, T>);
impl_bytes_conversion!(CrtGgsw<S, T>);
impl_zero!(CrtGgsw<S, T>);
impl_basic_operation_multiple_modulus!(CrtGgsw<S, T>);
impl_crt_ntt!(CrtGgsw<S, T>, DcrtGgsw);

impl<S, T> CrtGgsw<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_glev_mut(&mut self, crt_glev_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(crt_glev_len)
    }
}

impl<S, T> CrtGgsw<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_glev(&self, crt_glev_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(crt_glev_len)
    }
}
