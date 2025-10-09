use primus_integer::UnsignedInteger;

mod ops;

/// Unsigned integer modulus.
///
/// Just store the modulus value and only support some basic operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct UintModulus<T>(pub T);

impl<T: UnsignedInteger> UintModulus<T> {
    /// Creates a new [`UintModulus<T>`].
    #[inline(always)]
    pub fn new(value: T) -> Self {
        Self(value)
    }
}

impl<T: UnsignedInteger> primus_reduce::Modulus for UintModulus<T> {
    type ValueT = T;

    #[inline(always)]
    fn value(self) -> Option<Self::ValueT> {
        Some(self.0)
    }

    #[inline(always)]
    fn value_unchecked(self) -> Self::ValueT {
        self.0
    }

    #[inline(always)]
    fn minus_one(self) -> Self::ValueT {
        self.0 - T::ONE
    }
}
