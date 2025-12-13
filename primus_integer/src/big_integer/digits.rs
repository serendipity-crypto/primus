use std::slice::{Iter, IterMut};

pub trait Digits {
    type Limb: Copy;

    fn digits(&self) -> &[Self::Limb];

    /// Returns the number of elements.
    fn len(&self) -> usize;

    /// Returns an iterator.
    ///
    /// The iterator yields all items from start to end.
    fn iter<'a>(&'a self) -> Iter<'a, Self::Limb>;
}

pub trait DigitsMut: Digits {
    fn digits_mut(&mut self) -> &mut [<Self as Digits>::Limb];

    /// Returns an iterator that allows modifying each value.
    ///
    /// The iterator yields all items from start to end.
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, <Self as Digits>::Limb>;

    /// Fills `self` with elements by cloning `value`.
    fn fill(&mut self, value: <Self as Digits>::Limb);
}

impl<T: Copy> Digits for &[T] {
    type Limb = T;

    #[inline(always)]
    fn digits(&self) -> &[T] {
        self
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, Self::Limb> {
        <[T]>::iter(self)
    }
}

impl<T: Copy> Digits for &mut [T] {
    type Limb = T;

    #[inline(always)]
    fn digits(&self) -> &[T] {
        self
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, Self::Limb> {
        <[T]>::iter(self)
    }
}

impl<T: Copy> DigitsMut for &mut [T] {
    #[inline(always)]
    fn digits_mut(&mut self) -> &mut [T] {
        self
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, Self::Limb> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T) {
        <[T]>::fill(self, value);
    }
}

impl<T: Copy> Digits for Vec<T> {
    type Limb = T;

    #[inline(always)]
    fn digits(&self) -> &[T] {
        self.as_slice()
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <Vec<T>>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, Self::Limb> {
        <[T]>::iter(self)
    }
}

impl<T: Copy> DigitsMut for Vec<T> {
    #[inline(always)]
    fn digits_mut(&mut self) -> &mut [T] {
        self.as_mut_slice()
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, Self::Limb> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T) {
        <[T]>::fill(self, value);
    }
}

impl<T: Copy> Digits for Box<[T]> {
    type Limb = T;

    #[inline(always)]
    fn digits(&self) -> &[T] {
        self.as_ref()
    }

    #[inline(always)]
    fn len(&self) -> usize {
        <[T]>::len(self)
    }

    #[inline(always)]
    fn iter<'a>(&'a self) -> Iter<'a, Self::Limb> {
        <[T]>::iter(self)
    }
}

impl<T: Copy> DigitsMut for Box<[T]> {
    #[inline(always)]
    fn digits_mut(&mut self) -> &mut [T] {
        self.as_mut()
    }

    #[inline(always)]
    fn iter_mut<'a>(&'a mut self) -> IterMut<'a, Self::Limb> {
        <[T]>::iter_mut(self)
    }

    #[inline(always)]
    fn fill(&mut self, value: T) {
        <[T]>::fill(self, value);
    }
}
