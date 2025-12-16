use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use primus_integer::{Data, DataMut, RawData, UnsignedInteger};

use super::Polynomial;

impl<S, T, I: SliceIndex<[T]>> Index<I> for Polynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(self.as_ref(), index)
    }
}

impl<S, T, I: SliceIndex<[T]>> IndexMut<I> for Polynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(self.as_mut(), index)
    }
}

impl<S, T> AsRef<[T]> for Polynomial<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_slice()
    }
}

impl<S, T> AsMut<[T]> for Polynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut_slice()
    }
}
