use gcd::Xgcd;

use crate::{BorrowingSub, CarryingAdd, CarryingMul, DivRem, DivRemScalar, Integer, WideningMul};

/// An abstract over unsigned interger type.
pub trait UnsignedInteger:
    Integer
    + num_traits::Unsigned
    + CarryingAdd<CarryT = bool>
    + BorrowingSub<BorrowT = bool>
    + WideningMul
    + CarryingMul
    + DivRem
    + DivRemScalar
    + Xgcd
    + TryFrom<usize>
    + TryInto<usize>
{
    /// Returns `true` if and only if `self == 2^k` for some `k`.
    #[must_use]
    #[inline(always)]
    fn is_power_of_two(self) -> bool {
        self.count_ones() == 1
    }
}

macro_rules! impl_unsigned_integer {
    ($($t:ty)*) => ($(
        impl UnsignedInteger for $t {}
    )*)
}

impl_unsigned_integer! {u8 u16 u32 u64 u128 usize}
