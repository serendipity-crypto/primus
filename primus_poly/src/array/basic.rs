use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use num_traits::{ConstZero, Zero};
use primus_integer::{ByteCount, UnsignedInteger, size::Size};

use super::{ArrayBase, Data, DataMut, DataOwned, RawData};

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Returns an iterator that allows reading each value.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(<T as Zero>::is_zero)
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Returns an iterator that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.0.iter_mut()
    }

    #[inline]
    pub fn set_zero(&mut self) {
        self.0.fill(<T as ConstZero>::ZERO);
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    #[inline]
    pub fn zero(len: usize) -> Self {
        Self(S::new_array(<T as ConstZero>::ZERO, len))
    }
}

impl<S, T> Size for ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<S, T, I: SliceIndex<[T]>> Index<I> for ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data + Index<I, Output = I::Output>,
    T: UnsignedInteger,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        Index::index(&self.0, index)
    }
}

impl<S, T, I: SliceIndex<[T]>> IndexMut<I> for ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut + IndexMut<I, Output = I::Output>,
    T: UnsignedInteger,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut self.0, index)
    }
}

impl<S, T> AsRef<[T]> for ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<S, T> AsMut<[T]> for ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}

impl<'t, S, T> IntoIterator for &'t ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    type Item = &'t T;

    type IntoIter = core::slice::Iter<'t, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'t, S, T> IntoIterator for &'t mut ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    type Item = &'t mut T;

    type IntoIter = core::slice::IterMut<'t, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<S, T> IntoIterator for ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    type Item = T;

    type IntoIter = std::vec::IntoIter<T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        <S as DataOwned>::into_iter(self.0)
    }
}
