use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use super::CrtPolynomial;

use crate::Polynomial;

impl<T, I> IndexMut<I> for CrtPolynomial<T>
where
    I: SliceIndex<[Polynomial<T>]>,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut *self.polys, index)
    }
}

impl<T, I> Index<I> for CrtPolynomial<T>
where
    I: SliceIndex<[Polynomial<T>]>,
{
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&*self.polys, index)
    }
}

impl<T> AsRef<[Polynomial<T>]> for CrtPolynomial<T> {
    #[inline]
    fn as_ref(&self) -> &[Polynomial<T>] {
        self.polys.as_ref()
    }
}

impl<T> AsMut<[Polynomial<T>]> for CrtPolynomial<T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [Polynomial<T>] {
        self.polys.as_mut()
    }
}

impl<T> IntoIterator for CrtPolynomial<T> {
    type Item = Polynomial<T>;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a CrtPolynomial<T> {
    type Item = &'a Polynomial<T>;

    type IntoIter = core::slice::Iter<'a, Polynomial<T>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut CrtPolynomial<T> {
    type Item = &'a mut Polynomial<T>;

    type IntoIter = core::slice::IterMut<'a, Polynomial<T>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.iter_mut()
    }
}
