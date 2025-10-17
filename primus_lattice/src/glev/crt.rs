use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use super::DcrtGlev;

/// A representation of Module Learning with Errors (MLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct CrtGlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(CrtGlev<S, T>);
impl_bytes_conversion!(CrtGlev<S, T>);
impl_zero!(CrtGlev<S, T>);
impl_basic_operation_multiple_modulus!(CrtGlev<S, T>);
impl_crt_ntt!(CrtGlev<S, T>, DcrtGlev);

impl<S, T> CrtGlev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_glwe_mut(&mut self, crt_glwe_len: usize) -> std::slice::ChunksExactMut<'_, T> {
        self.data.chunks_exact_mut(crt_glwe_len)
    }
}

impl<S, T> CrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    pub fn iter_crt_glwe(&self, crt_glwe_len: usize) -> std::slice::ChunksExact<'_, T> {
        self.data.chunks_exact(crt_glwe_len)
    }
}
