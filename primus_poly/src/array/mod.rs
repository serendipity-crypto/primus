use std::ops::{Deref, DerefMut};

use primus_integer::UnsignedInteger;

mod basic;
mod random;

mod add;
mod mul;
mod neg;
mod sub;

pub type Array<T> = ArrayBase<Vec<T>>;

pub type ArrayRef<'a, T> = ArrayBase<&'a [T]>;

pub type ArrayMut<'a, T> = ArrayBase<&'a mut [T]>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArrayBase<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`ArrayBase<S, T>`].
    #[inline]
    pub fn new(data: S) -> Self {
        Self(data)
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    #[inline]
    pub fn from_slice(data: &[T]) -> Self {
        Self(S::from_slice(data))
    }

    #[inline(always)]
    pub fn from_vec(data: Vec<T>) -> Self {
        Self(S::from_vec(data))
    }

    #[inline]
    pub fn to_ref(&self) -> ArrayBase<&[T], T> {
        ArrayBase(self.0.to_ref())
    }
}

impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn copy_from_slice(&mut self, src: &[T]) {
        self.0.copy_from_slice(src);
    }

    #[inline]
    pub fn chunks_exact_mut(&mut self, chunk_size: usize) -> std::slice::ChunksExactMut<'_, T> {
        DataMut::chunks_exact_mut(&mut self.0, chunk_size)
    }
}
impl<S, T> ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Returns the number of elements.
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns `true` if `self` has a length of 0.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[inline]
    pub fn chunks_exact(&self, chunk_size: usize) -> std::slice::ChunksExact<'_, T> {
        Data::chunks_exact(&self.0, chunk_size)
    }
}

impl<S, T> FromIterator<T> for ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    #[inline]
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(<S as FromIterator<T>>::from_iter(iter))
    }
}

impl<S, T> Deref for ArrayBase<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    type Target = S;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<S, T> DerefMut for ArrayBase<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

pub trait RawData: Sized {
    type Elem: UnsignedInteger;
}

impl<T: UnsignedInteger> RawData for &[T] {
    type Elem = T;
}

impl<T: UnsignedInteger> RawData for &mut [T] {
    type Elem = T;
}

impl<T: UnsignedInteger> RawData for Vec<T> {
    type Elem = T;
}

pub trait Data: RawData + AsRef<[<Self as RawData>::Elem]> {
    /// Returns the number of elements.
    fn len(&self) -> usize;

    /// Returns `true` if `self` has a length of 0.
    fn is_empty(&self) -> bool;

    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end.
    fn iter<'a>(&'a self) -> core::slice::Iter<'a, <Self as RawData>::Elem>;

    fn chunks_exact<'a>(
        &'a self,
        chunk_size: usize,
    ) -> std::slice::ChunksExact<'a, <Self as RawData>::Elem>;

    /// Divides one slice into two at an index.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding
    /// the index `mid` itself) and the second will contain all
    /// indices from `[mid, len)` (excluding the index `len` itself).
    ///
    /// # Panics
    ///
    /// Panics if `mid > len`.
    fn split_at(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]);

    /// Divides one slice into two at an index, without doing bounds checking.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding the index `mid` itself)
    /// and the second will contain all indices from `[mid, len)` (excluding the index len itself).
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]);
}

impl<T: UnsignedInteger> Data for &[T] {
    #[inline]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        <[T]>::is_empty(self)
    }

    #[inline]
    fn iter<'a>(&'a self) -> core::slice::Iter<'a, <Self as RawData>::Elem> {
        <[T]>::iter(self)
    }

    #[inline]
    fn chunks_exact<'a>(
        &'a self,
        chunk_size: usize,
    ) -> std::slice::ChunksExact<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline]
    fn split_at(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        <[T]>::split_at(self, mid)
    }

    #[inline]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T: UnsignedInteger> Data for &mut [T] {
    #[inline]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        <[T]>::is_empty(self)
    }

    #[inline]
    fn iter<'a>(&'a self) -> core::slice::Iter<'a, <Self as RawData>::Elem> {
        <[T]>::iter(self)
    }

    #[inline]
    fn chunks_exact<'a>(
        &'a self,
        chunk_size: usize,
    ) -> std::slice::ChunksExact<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline]
    fn split_at(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        <[T]>::split_at(self, mid)
    }

    #[inline]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T: UnsignedInteger> Data for Vec<T> {
    #[inline]
    fn len(&self) -> usize {
        <Vec<T>>::len(self)
    }

    #[inline]
    fn is_empty(&self) -> bool {
        Vec::is_empty(self)
    }

    #[inline]
    fn iter<'a>(&'a self) -> core::slice::Iter<'a, <Self as RawData>::Elem> {
        <[T]>::iter(self)
    }

    #[inline]
    fn chunks_exact<'a>(
        &'a self,
        chunk_size: usize,
    ) -> std::slice::ChunksExact<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline]
    fn split_at(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        <[T]>::split_at(self, mid)
    }

    #[inline]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[Self::Elem], &[Self::Elem]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

pub trait DataMut: Data + AsMut<[<Self as RawData>::Elem]> {
    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end.
    fn iter_mut<'a>(&'a mut self) -> core::slice::IterMut<'a, <Self as RawData>::Elem>;

    fn chunks_exact_mut<'a>(
        &'a mut self,
        chunk_size: usize,
    ) -> std::slice::ChunksExactMut<'a, <Self as RawData>::Elem>;

    /// Fills `self` with elements by cloning `value`.
    fn fill(&mut self, value: <Self as RawData>::Elem);

    /// Copies all elements from `src` into `self`, using a memcpy.
    ///
    /// The length of `src` must be the same as `self`.
    fn copy_from_slice(&mut self, src: &[Self::Elem]);

    fn split_at_mut(&mut self, mid: usize) -> (&mut [Self::Elem], &mut [Self::Elem]);

    /// Divides one mutable slice into two at an index, without doing bounds checking.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding the index `mid` itself)
    /// and the second will contain all indices from `[mid, len)` (excluding the index len itself).
    unsafe fn split_at_mut_unchecked(
        &mut self,
        mid: usize,
    ) -> (&mut [Self::Elem], &mut [Self::Elem]);
}

impl<T: UnsignedInteger> DataMut for &mut [T] {
    #[inline]
    fn iter_mut<'a>(&'a mut self) -> core::slice::IterMut<'a, <Self as RawData>::Elem> {
        <[T]>::iter_mut(self)
    }

    #[inline]
    fn fill(&mut self, value: T) {
        <[T]>::fill(self, value);
    }

    #[inline]
    fn copy_from_slice(&mut self, src: &[Self::Elem]) {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline]
    fn chunks_exact_mut<'a>(
        &'a mut self,
        chunk_size: usize,
    ) -> std::slice::ChunksExactMut<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [Self::Elem], &mut [Self::Elem]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(
        &mut self,
        mid: usize,
    ) -> (&mut [Self::Elem], &mut [Self::Elem]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

impl<T: UnsignedInteger> DataMut for Vec<T> {
    #[inline]
    fn iter_mut<'a>(&'a mut self) -> core::slice::IterMut<'a, <Self as RawData>::Elem> {
        <[T]>::iter_mut(self)
    }

    #[inline]
    fn fill(&mut self, value: T) {
        <[T]>::fill(self, value);
    }

    #[inline]
    fn copy_from_slice(&mut self, src: &[Self::Elem]) {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline]
    fn chunks_exact_mut<'a>(
        &'a mut self,
        chunk_size: usize,
    ) -> std::slice::ChunksExactMut<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [Self::Elem], &mut [Self::Elem]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(
        &mut self,
        mid: usize,
    ) -> (&mut [Self::Elem], &mut [Self::Elem]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

pub trait DataOwned: DataMut + FromIterator<Self::Elem> {
    /// Creates array.
    fn new_array(value: Self::Elem, len: usize) -> Self;

    fn from_slice(data: &[Self::Elem]) -> Self;

    fn from_vec(data: Vec<Self::Elem>) -> Self;

    /// Creates a consuming iterator, that is, one that moves each value out of the vector (from start to end).
    fn into_iter(self) -> std::vec::IntoIter<<Self as RawData>::Elem>;

    fn to_ref(&self) -> &[Self::Elem];
}

impl<T: UnsignedInteger> DataOwned for Vec<T> {
    #[inline(always)]
    fn new_array(value: Self::Elem, len: usize) -> Self {
        vec![value; len]
    }

    #[inline]
    fn from_slice(data: &[Self::Elem]) -> Self {
        data.to_vec()
    }

    #[inline(always)]
    fn from_vec(data: Vec<Self::Elem>) -> Self {
        data
    }

    #[inline]
    fn into_iter(self) -> std::vec::IntoIter<<Self as RawData>::Elem> {
        <Vec<T> as IntoIterator>::into_iter(self)
    }

    #[inline]
    fn to_ref(&self) -> &[Self::Elem] {
        self
    }
}
