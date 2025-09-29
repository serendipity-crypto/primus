use core::{
    fmt::{Debug, Display},
    ops::{
        BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Not, Shl, ShlAssign, Shr,
        ShrAssign,
    },
};

use num_traits::{ConstOne, ConstZero, FromBytes, MulAdd, MulAddAssign, NumAssign, Pow, ToBytes};
use primus_utils::ByteCount;
use rand::distr::uniform::SampleUniform;
use serde::{Deserialize, Serialize};
use zeroize::Zeroize;

use crate::integer_traits::*;

/// An abstract over interger type.
pub trait Integer:
    Sized
    + Send
    + Sync
    + Clone
    + Copy
    + Default
    + PartialOrd
    + Ord
    + PartialEq
    + Eq
    + Debug
    + Display
    + Bits
    + ByteCount
    + ToBytes
    + FromBytes
    + ConstZero
    + ConstOne
    + ConstTwo
    + ConstBounded
    + AsCast
    + AsFrom<bool>
    + NumAssign
    + WrappingAdd
    + WrappingSub
    + WrappingNeg
    + WrappingMul
    + WrappingShl
    + WrappingShr
    + OverflowingAdd
    + OverflowingSub
    + OverflowingMul
    + CheckedAdd
    + CheckedSub
    + CheckedMul
    + CheckedDiv
    + CheckedNeg
    + CheckedRem
    + CheckedShl
    + CheckedShr
    + MulAdd
    + MulAddAssign
    + Not<Output = Self>
    + BitAnd<Output = Self>
    + BitOr<Output = Self>
    + BitXor<Output = Self>
    + BitAndAssign
    + BitOrAssign
    + BitXorAssign
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
    + ShlAssign<u32>
    + ShrAssign<u32>
    + Pow<u32, Output = Self>
    + Pow<usize, Output = Self>
    + SampleUniform<Sampler: Copy>
    + Serialize
    + for<'de> Deserialize<'de>
    + Zeroize
{
}

macro_rules! empty_trait_impl {
    ($name:ident for $($t:ty)*) => ($(
        impl $name for $t {}
    )*)
}

// `isize` is not supported by `rand::distr::uniform::SampleUniform`.
empty_trait_impl!(Integer for u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128);
