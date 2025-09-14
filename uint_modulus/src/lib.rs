use integer::UnsignedInteger;

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
