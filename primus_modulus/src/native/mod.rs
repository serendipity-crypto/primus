use core::{marker::PhantomData, panic};

use primus_integer::UnsignedInteger;

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

impl<T: UnsignedInteger> primus_reduce::Modulus for NativeModulus<T> {
    type ValueT = T;

    #[inline(always)]
    fn value(self) -> Option<Self::ValueT> {
        None
    }

    #[inline(always)]
    fn value_unchecked(self) -> Self::ValueT {
        panic!("The value of theNative Modulus can not be represented.");
    }

    #[inline(always)]
    fn minus_one(self) -> Self::ValueT {
        T::MAX
    }
}
