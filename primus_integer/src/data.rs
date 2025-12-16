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
    fn is_empty(&self) -> bool;

    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end.
    fn iter<'a>(&'a self) -> Iter<'a, Self::Elem>;
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
    fn is_empty(&self) -> bool {
        <[T]>::is_empty(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
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
    fn is_empty(&self) -> bool {
        <[T]>::is_empty(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
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
    fn is_empty(&self) -> bool {
        Vec::is_empty(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
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
    fn is_empty(&self) -> bool {
        <[T]>::is_empty(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, T> {
        <[T]>::iter(self)
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
}
