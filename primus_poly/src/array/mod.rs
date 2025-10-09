mod basic;
mod random;

mod add;
mod mul;
mod neg;
mod sub;

#[derive(Debug, Clone)]
pub struct Array<T>(pub Vec<T>);

impl<T> Array<T> {
    #[inline]
    pub fn new(data: Vec<T>) -> Self {
        Self(data)
    }

    #[inline]
    pub fn to_ref<'a>(&'a self) -> ArrayRef<'a, T> {
        ArrayRef(&self.0)
    }

    #[inline]
    pub fn to_mut<'a>(&'a mut self) -> ArrayMut<'a, T> {
        ArrayMut(&mut self.0)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ArrayRef<'a, T>(pub &'a [T]);

impl<'a, T> ArrayRef<'a, T> {
    #[inline]
    pub fn new(data: &'a [T]) -> Self {
        Self(data)
    }
}

#[derive(Debug)]
pub struct ArrayMut<'a, T>(pub &'a mut [T]);

impl<'a, T> ArrayMut<'a, T> {
    #[inline]
    pub fn new(data: &'a mut [T]) -> Self {
        Self(data)
    }

    #[inline]
    pub fn to_ref(self) -> ArrayRef<'a, T> {
        ArrayRef(&*self.0)
    }
}
