#![cfg_attr(all(feature = "nightly", feature = "simd"), feature(portable_simd))]

use core::fmt::Display;

pub use primus_integer as integer;
pub use primus_reduce as reduce;
pub use primus_uint_modulus as uint_modulus;

use integer::{DivRemScalar, UnsignedInteger};

#[cfg(feature = "derive")]
pub use primus_barrett_derive::Barrett;

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
    pub fn new(value: T) -> Self {
        assert!(value > T::ONE, "modulus can't be 0 or 1.");

        let bit_count = T::BITS - value.leading_zeros();
        assert!(bit_count < T::BITS - 1, "modulus is too large.");

        let mut quotient = [T::ZERO; 3];
        let _rem = DivRemScalar::div_rem_scalar(&[T::ZERO, T::ZERO, T::ONE], value, &mut quotient);

        Self {
            value,
            ratio: [quotient[0], quotient[1]],
        }
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

impl<T: UnsignedInteger> reduce::Modulus for BarrettModulus<T> {
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
    use reduce::FieldContext;

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
