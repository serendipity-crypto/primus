use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use super::DcrtPolynomial;

use crate::NttPolynomial;

impl<T, I> IndexMut<I> for DcrtPolynomial<T>
where
    I: SliceIndex<[NttPolynomial<T>]>,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut *self.polys, index)
    }
}

impl<T, I> Index<I> for DcrtPolynomial<T>
where
    I: SliceIndex<[NttPolynomial<T>]>,
{
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&*self.polys, index)
    }
}

impl<T> AsRef<[NttPolynomial<T>]> for DcrtPolynomial<T> {
    #[inline]
    fn as_ref(&self) -> &[NttPolynomial<T>] {
        self.polys.as_ref()
    }
}

impl<T> AsMut<[NttPolynomial<T>]> for DcrtPolynomial<T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [NttPolynomial<T>] {
        self.polys.as_mut()
    }
}

impl<T> IntoIterator for DcrtPolynomial<T> {
    type Item = NttPolynomial<T>;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a DcrtPolynomial<T> {
    type Item = &'a NttPolynomial<T>;

    type IntoIter = core::slice::Iter<'a, NttPolynomial<T>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut DcrtPolynomial<T> {
    type Item = &'a mut NttPolynomial<T>;

    type IntoIter = core::slice::IterMut<'a, NttPolynomial<T>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.polys.iter_mut()
    }
}
