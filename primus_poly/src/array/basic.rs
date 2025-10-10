use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use num_traits::{ConstZero, Zero};
use primus_integer::{ByteCount, UnsignedInteger, size::Size};
use primus_reduce::ops::ReduceMulAdd;

use super::{ArrayBase, Data, DataMut, DataOwned, RawData};

impl<S, T> ArrayBase<S>
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

    /// Evaluate the polynomial with the value `x`.
    #[inline]
    pub fn evaluate<M>(&self, x: T, modulus: M) -> T
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.0.iter().rev().fold(<T as ConstZero>::ZERO, |acc, &a| {
            modulus.reduce_mul_add(acc, x, a)
        })
    }
}

impl<S, T> ArrayBase<S>
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

impl<S, T> ArrayBase<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    #[inline]
    pub fn zero(len: usize) -> Self {
        Self(S::new_array(<T as ConstZero>::ZERO, len))
    }
}

impl<S, T> Size for ArrayBase<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<S, T, I: SliceIndex<[T]>> Index<I> for ArrayBase<S>
where
    S: RawData<Elem = T> + Data + Index<I, Output = I::Output>,
    T: UnsignedInteger,
{
    type Output = I::Output;

    fn index(&self, index: I) -> &Self::Output {
        Index::index(&self.0, index)
    }
}

impl<S, T, I: SliceIndex<[T]>> IndexMut<I> for ArrayBase<S>
where
    S: RawData<Elem = T> + DataMut + IndexMut<I, Output = I::Output>,
    T: UnsignedInteger,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut self.0, index)
    }
}

impl<S, T> AsRef<[T]> for ArrayBase<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<S, A> AsMut<[A]> for ArrayBase<S>
where
    S: RawData<Elem = A> + DataMut,
    A: UnsignedInteger,
{
    #[inline]
    fn as_mut(&mut self) -> &mut [A] {
        self.0.as_mut()
    }
}

impl<'t, S, A> IntoIterator for &'t ArrayBase<S>
where
    S: RawData<Elem = A> + Data,
    A: UnsignedInteger,
{
    type Item = &'t A;

    type IntoIter = core::slice::Iter<'t, A>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'t, S, A> IntoIterator for &'t mut ArrayBase<S>
where
    S: RawData<Elem = A> + DataMut,
    A: UnsignedInteger,
{
    type Item = &'t mut A;

    type IntoIter = core::slice::IterMut<'t, A>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<S, A> IntoIterator for ArrayBase<S>
where
    S: RawData<Elem = A> + DataOwned,
    A: UnsignedInteger,
{
    type Item = A;

    type IntoIter = std::vec::IntoIter<A>;

    fn into_iter(self) -> Self::IntoIter {
        <S as DataOwned>::into_iter(self.0)
    }
}
