use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use primus_integer::UnsignedInteger;

use crate::{Data, DataMut, RawData};

use super::Polynomial;

impl<S, T, I: SliceIndex<[T]>> Index<I> for Polynomial<S, T>
where
    S: RawData<Elem = T> + Data + Index<I, Output = I::Output>,
    T: UnsignedInteger,
{
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&self.0, index)
    }
}

impl<S, T, I: SliceIndex<[T]>> IndexMut<I> for Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut + IndexMut<I, Output = I::Output>,
    T: UnsignedInteger,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut self.0, index)
    }
}

impl<S, T> AsRef<[T]> for Polynomial<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<S, T> AsMut<[T]> for Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}
