use core::fmt::Display;

use crate::integer::{DivRemScalar, UnsignedInteger};

mod ops;
#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

#[cfg(all(feature = "nightly", feature = "simd"))]
pub use simd::SimdBarrettModulus;

/// A modulus, using barrett reduction algorithm.
///
/// The struct stores the modulus number and some precomputed
/// data. Here, `b` = 2^T::BITS
///
/// It's efficient if many reductions are performed with a single modulus.
#[derive(Debug, Clone, Copy)]
pub struct BarrettModulus<T: UnsignedInteger> {
    /// the value to indicate the modulus
    value: T,
    /// ratio `µ` = floor(b²/value)
    ratio: [T; 2],
}

impl<T: UnsignedInteger> BarrettModulus<T> {
    /// Creates a new [`BarrettModulus<T>`] with the given value.
    ///
    /// # Panics
    ///
    /// Panics if `value ≤ 1` or if the bit-width of `value` is too large
    /// (≥ `T::BITS - 1`).  For a fallible variant see [`try_new`](Self::try_new).
    pub fn new(value: T) -> Self {
        assert!(value > T::ONE, "modulus can't be 0 or 1.");
        let bit_count = T::BITS - value.leading_zeros();
        assert!(bit_count < T::BITS - 1, "modulus is too large.");
        Self::new_unchecked(value)
    }

    /// Creates a new [`BarrettModulus<T>`] without validating the modulus value.
    ///
    /// # Correctness
    ///
    /// `value` must satisfy `1 < value < 2^(T::BITS - 1)`.
    #[inline]
    pub fn new_unchecked(value: T) -> Self {
        let mut quotient = [T::ZERO; 3];
        let _rem = DivRemScalar::div_rem_scalar(&[T::ZERO, T::ZERO, T::ONE], value, &mut quotient);
        Self {
            value,
            ratio: [quotient[0], quotient[1]],
        }
    }

    /// Fallible constructor returning `None` if the value is out of range.
    #[inline]
    pub fn try_new(value: T) -> Option<Self> {
        if value <= T::ONE {
            return None;
        }
        let bit_count = T::BITS - value.leading_zeros();
        if bit_count >= T::BITS - 1 {
            return None;
        }
        Some(Self::new_unchecked(value))
    }

    /// Returns the value of this [`BarrettModulus<T>`].
    #[inline]
    pub const fn value(&self) -> T {
        self.value
    }

    /// Returns the ratio of this [`BarrettModulus<T>`].
    #[inline]
    pub const fn ratio(&self) -> [T; 2] {
        self.ratio
    }
}

impl<T: UnsignedInteger> Display for BarrettModulus<T> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<T: UnsignedInteger> primus_reduce::Modulus for BarrettModulus<T> {
    type ValueT = T;

    #[inline]
    fn value(self) -> Option<Self::ValueT> {
        Some(self.value)
    }

    #[inline(always)]
    fn value_unchecked(self) -> Self::ValueT {
        self.value
    }

    #[inline(always)]
    fn minus_one(self) -> Self::ValueT {
        self.value - T::ONE
    }
}

#[cfg(test)]
mod tests {
    use primus_reduce::FieldContext;

    use super::*;

    fn field_trait<T: UnsignedInteger, M: FieldContext<T>>(_modulus: M) {}

    #[test]
    fn test_trait() {
        field_trait(<BarrettModulus<u8>>::new(61));
        field_trait(<BarrettModulus<u16>>::new(12289));
        field_trait(<BarrettModulus<u32>>::new(536813569));
        field_trait(<BarrettModulus<u64>>::new(4611686018427322369));
    }
}
