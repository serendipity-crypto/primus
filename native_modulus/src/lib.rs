use core::marker::PhantomData;

use integer::UnsignedInteger;

mod ops;

/// Natvie modulus.
///
/// - For `u8`, this type acts as `2⁸`
/// - For `u16`, this type acts as `2¹⁶`
/// - For `u32`, this type acts as `2³²`
/// - For `u64`, this type acts as `2⁶⁴`
/// - For `u128`, this type acts as `2¹²⁸`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct NativeModulus<T: UnsignedInteger> {
    phantom: PhantomData<T>,
}

impl<T: UnsignedInteger> Default for NativeModulus<T> {
    #[inline(always)]
    fn default() -> Self {
        Self::new()
    }
}

impl<T: UnsignedInteger> NativeModulus<T> {
    /// Creates a new [`NativeModulus<T>`].
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            phantom: PhantomData,
        }
    }
}
