use crate::ByteCount;

/// A trait for the size of a value.
pub trait Size {
    /// Returns the size of the pointed-to value in bytes.
    fn bytes_count(&self) -> usize;
}

impl<T: ByteCount> Size for Vec<T> {
    #[inline]
    fn bytes_count(&self) -> usize {
        self.len() * T::BYTES_COUNT
    }
}

impl<T: ByteCount> Size for &[T] {
    #[inline]
    fn bytes_count(&self) -> usize {
        self.len() * T::BYTES_COUNT
    }
}
