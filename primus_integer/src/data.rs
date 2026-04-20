use std::slice::{Iter, IterMut};

pub trait RawData: Sized {
    type Elem;
}

pub trait Data: RawData {
    /// Returns a slice containing the entire data.
    fn as_slice(&self) -> &[Self::Elem];

    /// Returns the number of elements.
    fn len(&self) -> usize;

    /// Returns `true` if `self` has a length of 0.
    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end.
    fn iter<'a>(&'a self) -> Iter<'a, Self::Elem>;

    /// Returns an iterator over `chunk_size` elements of the slice at a time, starting at the
    /// beginning of the slice.
    ///
    /// The chunks are slices and do not overlap. If `chunk_size` does not divide the length of the
    /// slice, then the last up to `chunk_size-1` elements will be omitted and can be retrieved
    /// from the `remainder` function of the iterator.
    ///
    /// Due to each chunk having exactly `chunk_size` elements, the compiler can often optimize the
    /// resulting code better than in the case of [`chunks`].
    ///
    /// # Panics
    ///
    /// Panics if `chunk_size` is zero.
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
    fn split_at(&self, mid: usize) -> (&[<Self as RawData>::Elem], &[<Self as RawData>::Elem]);

    /// Divides one slice into two at an index, without doing bounds checking.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding
    /// the index `mid` itself) and the second will contain all
    /// indices from `[mid, len)` (excluding the index `len` itself).
    ///
    /// For a safe alternative see [`split_at`].
    ///
    /// # Safety
    ///
    /// Calling this method with an out-of-bounds index is *[undefined behavior]*
    /// even if the resulting reference is not used. The caller has to ensure that
    /// `0 <= mid <= self.len()`.
    ///
    /// [undefined behavior]: https://doc.rust-lang.org/reference/behavior-considered-undefined.html
    unsafe fn split_at_unchecked(
        &self,
        mid: usize,
    ) -> (&[<Self as RawData>::Elem], &[<Self as RawData>::Elem]);
}

pub trait DataMut: Data {
    /// Returns a mutable slice containing the entire data.
    fn as_mut_slice(&mut self) -> &mut [<Self as RawData>::Elem];

    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end.
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, <Self as RawData>::Elem>;

    /// Fills `self` with elements by cloning `value`.
    fn fill(&mut self, value: <Self as RawData>::Elem)
    where
        <Self as RawData>::Elem: Clone;

    /// Copies all elements from `src` into `self`, using a memcpy.
    ///
    /// The length of `src` must be the same as `self`.
    fn copy_from_slice(&mut self, src: &[<Self as RawData>::Elem])
    where
        <Self as RawData>::Elem: Copy;

    fn chunks_exact_mut<'a>(
        &'a mut self,
        chunk_size: usize,
    ) -> std::slice::ChunksExactMut<'a, <Self as RawData>::Elem>;

    /// Divides one mutable slice into two at an index.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding
    /// the index `mid` itself) and the second will contain all
    /// indices from `[mid, len)` (excluding the index `len` itself).
    ///
    /// # Panics
    ///
    /// Panics if `mid > len`.
    fn split_at_mut(
        &mut self,
        mid: usize,
    ) -> (
        &mut [<Self as RawData>::Elem],
        &mut [<Self as RawData>::Elem],
    );

    /// Divides one mutable slice into two at an index, without doing bounds checking.
    ///
    /// The first will contain all indices from `[0, mid)` (excluding
    /// the index `mid` itself) and the second will contain all
    /// indices from `[mid, len)` (excluding the index `len` itself).
    ///
    /// For a safe alternative see [`split_at_mut`].
    ///
    /// # Safety
    ///
    /// Calling this method with an out-of-bounds index is *[undefined behavior]*
    /// even if the resulting reference is not used. The caller has to ensure that
    /// `0 <= mid <= self.len()`.
    ///
    /// [`split_at_mut`]: std::slice::split_at_mut
    /// [undefined behavior]: https://doc.rust-lang.org/reference/behavior-considered-undefined.html
    unsafe fn split_at_mut_unchecked(
        &mut self,
        mid: usize,
    ) -> (
        &mut [<Self as RawData>::Elem],
        &mut [<Self as RawData>::Elem],
    );
}

pub trait DataOwned: Data + FromIterator<<Self as RawData>::Elem> {
    type IntoIter: Iterator<Item = <Self as RawData>::Elem>;

    fn from_slice(data: &[<Self as RawData>::Elem]) -> Self
    where
        <Self as RawData>::Elem: Clone;

    fn from_vec(data: Vec<<Self as RawData>::Elem>) -> Self;

    /// Creates a consuming iterator, that is, one that moves each value out of the vector (from start to end).
    fn into_iter(self) -> Self::IntoIter;
}

impl<T> RawData for &[T] {
    type Elem = T;
}

impl<T> Data for &[T] {
    #[inline(always)]
    fn as_slice(&self) -> &[T] {
        self
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
    }

    #[inline(always)]
    fn chunks_exact<'a>(&'a self, chunk_size: usize) -> std::slice::ChunksExact<'a, T> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline(always)]
    fn split_at(&self, mid: usize) -> (&[T], &[T]) {
        <[T]>::split_at(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T> RawData for &mut [T] {
    type Elem = T;
}

impl<T> Data for &mut [T] {
    #[inline(always)]
    fn as_slice(&self) -> &[T] {
        self
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
    }

    #[inline(always)]
    fn chunks_exact<'a>(&'a self, chunk_size: usize) -> std::slice::ChunksExact<'a, T> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline(always)]
    fn split_at(&self, mid: usize) -> (&[T], &[T]) {
        <[T]>::split_at(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T> DataMut for &mut [T] {
    #[inline(always)]
    fn as_mut_slice(&mut self) -> &mut [T] {
        self
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        <[T]>::fill(self, value);
    }

    #[inline(always)]
    fn copy_from_slice(&mut self, src: &[T])
    where
        T: Copy,
    {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline(always)]
    fn chunks_exact_mut<'a>(&'a mut self, chunk_size: usize) -> std::slice::ChunksExactMut<'a, T> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline(always)]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

impl<T, const N: usize> RawData for [T; N] {
    type Elem = T;
}

impl<T, const N: usize> Data for [T; N] {
    #[inline(always)]
    fn as_slice(&self) -> &[Self::Elem] {
        self
    }

    #[inline(always)]
    fn len(&self) -> usize {
        N
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, Self::Elem> {
        <[T]>::iter(self)
    }

    #[inline(always)]
    fn chunks_exact<'a>(
        &'a self,
        chunk_size: usize,
    ) -> std::slice::ChunksExact<'a, <Self as RawData>::Elem> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline(always)]
    fn split_at(&self, mid: usize) -> (&[<Self as RawData>::Elem], &[<Self as RawData>::Elem]) {
        <[T]>::split_at(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_unchecked(
        &self,
        mid: usize,
    ) -> (&[<Self as RawData>::Elem], &[<Self as RawData>::Elem]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T, const N: usize> DataMut for [T; N] {
    #[inline(always)]
    fn as_mut_slice(&mut self) -> &mut [T] {
        self
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        <[T]>::fill(self, value);
    }

    #[inline(always)]
    fn copy_from_slice(&mut self, src: &[T])
    where
        T: Copy,
    {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline(always)]
    fn chunks_exact_mut<'a>(&'a mut self, chunk_size: usize) -> std::slice::ChunksExactMut<'a, T> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline(always)]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

impl<T> RawData for Vec<T> {
    type Elem = T;
}

impl<T> Data for Vec<T> {
    #[inline(always)]
    fn as_slice(&self) -> &[T] {
        self.as_slice()
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <Vec<T>>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
    }

    #[inline(always)]
    fn chunks_exact<'a>(&'a self, chunk_size: usize) -> std::slice::ChunksExact<'a, T> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline(always)]
    fn split_at(&self, mid: usize) -> (&[T], &[T]) {
        <[T]>::split_at(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T> DataMut for Vec<T> {
    #[inline(always)]
    fn as_mut_slice(&mut self) -> &mut [T] {
        self.as_mut_slice()
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        <[T]>::fill(self, value);
    }

    #[inline(always)]
    fn copy_from_slice(&mut self, src: &[T])
    where
        T: Copy,
    {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline(always)]
    fn chunks_exact_mut<'a>(&'a mut self, chunk_size: usize) -> std::slice::ChunksExactMut<'a, T> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline(always)]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

impl<T> DataOwned for Vec<T> {
    type IntoIter = std::vec::IntoIter<T>;

    #[inline(always)]
    fn from_slice(data: &[T]) -> Self
    where
        T: Clone,
    {
        data.to_vec()
    }

    #[inline(always)]
    fn from_vec(data: Vec<T>) -> Self {
        data
    }

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        <Vec<T> as IntoIterator>::into_iter(self)
    }
}

impl<T> RawData for Box<[T]> {
    type Elem = T;
}

impl<T> Data for Box<[T]> {
    #[inline(always)]
    fn as_slice(&self) -> &[T] {
        self.as_ref()
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
    }

    #[inline(always)]
    fn chunks_exact<'a>(&'a self, chunk_size: usize) -> std::slice::ChunksExact<'a, T> {
        <[T]>::chunks_exact(self, chunk_size)
    }

    #[inline(always)]
    fn split_at(&self, mid: usize) -> (&[T], &[T]) {
        <[T]>::split_at(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_unchecked(&self, mid: usize) -> (&[T], &[T]) {
        unsafe { <[T]>::split_at_unchecked(self, mid) }
    }
}

impl<T> DataMut for Box<[T]> {
    #[inline(always)]
    fn as_mut_slice(&mut self) -> &mut [T] {
        self.as_mut()
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, T> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T)
    where
        T: Clone,
    {
        <[T]>::fill(self, value);
    }

    #[inline(always)]
    fn copy_from_slice(&mut self, src: &[T])
    where
        T: Copy,
    {
        <[T]>::copy_from_slice(self, src);
    }

    #[inline(always)]
    fn chunks_exact_mut<'a>(&'a mut self, chunk_size: usize) -> std::slice::ChunksExactMut<'a, T> {
        <[T]>::chunks_exact_mut(self, chunk_size)
    }

    #[inline(always)]
    fn split_at_mut(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        <[T]>::split_at_mut(self, mid)
    }

    #[inline(always)]
    unsafe fn split_at_mut_unchecked(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        unsafe { <[T]>::split_at_mut_unchecked(self, mid) }
    }
}

impl<T> DataOwned for Box<[T]> {
    type IntoIter = <Box<[T]> as IntoIterator>::IntoIter;

    #[inline(always)]
    fn from_slice(data: &[T]) -> Self
    where
        T: Clone,
    {
        data.to_vec().into_boxed_slice()
    }

    #[inline(always)]
    fn from_vec(data: Vec<T>) -> Self {
        data.into_boxed_slice()
    }

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        <Box<[T]> as IntoIterator>::into_iter(self)
    }
}
