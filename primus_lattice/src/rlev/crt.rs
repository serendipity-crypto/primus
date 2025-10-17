use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::DcrtRlev;

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct CrtRlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtRlev<S, T>);
impl_bytes_conversion!(CrtRlev<S, T>);
impl_zero!(CrtRlev<S, T>);
impl_basic_operation_multiple_modulus!(CrtRlev<S, T>);
impl_crt_ntt!(CrtRlev<S, T>, DcrtRlev);

impl<S, T> CrtRlev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_rlwe_mut(&mut self, crt_rlwe_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(crt_rlwe_len)
    }
}

impl<S, T> CrtRlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_rlwe(&self, crt_rlwe_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(crt_rlwe_len)
    }
}
