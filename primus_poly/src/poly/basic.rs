use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use super::{Polynomial, PolynomialRef, PolynomialRefMut};

impl<T, I: SliceIndex<[T]>> Index<I> for Polynomial<T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&*self.data, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> Index<I> for PolynomialRef<'p, T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(self.data, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> Index<I> for PolynomialRefMut<'p, T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(self.data, index)
    }
}

impl<T, I: SliceIndex<[T]>> IndexMut<I> for Polynomial<T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut *self.data, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> IndexMut<I> for PolynomialRefMut<'p, T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(self.data, index)
    }
}

impl<T> AsRef<[T]> for Polynomial<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<'p, T> AsRef<[T]> for PolynomialRef<'p, T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.data
    }
}

impl<'p, T> AsRef<[T]> for PolynomialRefMut<'p, T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.data
    }
}

impl<T> AsMut<[T]> for Polynomial<T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.data.as_mut()
    }
}

impl<'p, T> AsMut<[T]> for PolynomialRefMut<'p, T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.data
    }
}

impl<T> IntoIterator for Polynomial<T> {
    type Item = T;

    type IntoIter = std::vec::IntoIter<T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Polynomial<T> {
    type Item = &'a T;

    type IntoIter = core::slice::Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut Polynomial<T> {
    type Item = &'a mut T;

    type IntoIter = core::slice::IterMut<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.data.iter_mut()
    }
}

impl<'p, T> IntoIterator for PolynomialRef<'p, T> {
    type Item = &'p T;

    type IntoIter = core::slice::Iter<'p, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter()
    }
}

impl<'p, T> IntoIterator for PolynomialRefMut<'p, T> {
    type Item = &'p mut T;

    type IntoIter = core::slice::IterMut<'p, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.data.iter_mut()
    }
}
