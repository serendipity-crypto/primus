use core::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use num_traits::{ConstZero, Zero};
use primus_integer::{ByteCount, UnsignedInteger, size::Size};
use primus_reduce::ops::ReduceMulAdd;

use super::{Array, ArrayMut, ArrayRef};

impl<T> Array<T> {
    /// Returns an iterator that allows reading each value.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }

    /// Returns an iterator that allows modifying each value.
    #[inline]
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.0.iter_mut()
    }
}

impl<T: Copy + ConstZero> Array<T> {
    #[inline]
    pub fn zero(len: usize) -> Self {
        Self(vec![<T as ConstZero>::ZERO; len])
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(<T as Zero>::is_zero)
    }

    #[inline]
    pub fn set_zero(&mut self) {
        self.0.fill(<T as ConstZero>::ZERO);
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

impl<'a, T> ArrayRef<'a, T> {
    /// Returns an iterator that allows reading each value.
    #[inline]
    pub fn iter(&self) -> std::slice::Iter<'a, T> {
        self.0.iter()
    }
}

impl<'a, T: Copy + ConstZero> ArrayRef<'a, T> {
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

impl<'a, T> ArrayMut<'a, T> {
    /// Returns an iterator that allows reading each value.
    #[inline]
    pub fn iter<'b>(&'b self) -> std::slice::Iter<'b, T> {
        self.0.iter()
    }

    /// Returns an iterator that allows modifying each value.
    #[inline]
    pub fn iter_mut<'b>(&'b mut self) -> std::slice::IterMut<'b, T> {
        self.0.iter_mut()
    }
}

impl<'a, T: Copy + ConstZero> ArrayMut<'a, T> {
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(<T as Zero>::is_zero)
    }

    #[inline]
    pub fn set_zero(&mut self) {
        self.0.fill(<T as ConstZero>::ZERO);
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

impl<T: UnsignedInteger> Size for Array<T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<'a, T: UnsignedInteger> Size for ArrayRef<'a, T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<'a, T: UnsignedInteger> Size for ArrayMut<'a, T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T, I: SliceIndex<[T]>> Index<I> for Array<T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&*self.0, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> Index<I> for ArrayRef<'p, T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(self.0, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> Index<I> for ArrayMut<'p, T> {
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(self.0, index)
    }
}

impl<T, I: SliceIndex<[T]>> IndexMut<I> for Array<T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut *self.0, index)
    }
}

impl<'p, T, I: SliceIndex<[T]>> IndexMut<I> for ArrayMut<'p, T> {
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(self.0, index)
    }
}

impl<T> AsRef<[T]> for Array<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0.as_ref()
    }
}

impl<'p, T> AsRef<[T]> for ArrayRef<'p, T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0
    }
}

impl<'p, T> AsRef<[T]> for ArrayMut<'p, T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        self.0
    }
}

impl<T> AsMut<[T]> for Array<T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0.as_mut()
    }
}

impl<'p, T> AsMut<[T]> for ArrayMut<'p, T> {
    #[inline]
    fn as_mut(&mut self) -> &mut [T] {
        self.0
    }
}

impl<T> IntoIterator for Array<T> {
    type Item = T;

    type IntoIter = std::vec::IntoIter<T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a, T> IntoIterator for &'a Array<T> {
    type Item = &'a T;

    type IntoIter = core::slice::Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut Array<T> {
    type Item = &'a mut T;

    type IntoIter = core::slice::IterMut<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'a, T> IntoIterator for ArrayRef<'a, T> {
    type Item = &'a T;

    type IntoIter = core::slice::Iter<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'b, 'a: 'b, T> IntoIterator for &'b ArrayRef<'a, T> {
    type Item = &'b T;

    type IntoIter = core::slice::Iter<'b, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a, T> IntoIterator for ArrayMut<'a, T> {
    type Item = &'a mut T;

    type IntoIter = core::slice::IterMut<'a, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<'b, 'a: 'b, T> IntoIterator for &'b ArrayMut<'a, T> {
    type Item = &'b T;

    type IntoIter = core::slice::Iter<'b, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'b, 'a: 'b, T> IntoIterator for &'b mut ArrayMut<'a, T> {
    type Item = &'b mut T;

    type IntoIter = core::slice::IterMut<'b, T>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}
