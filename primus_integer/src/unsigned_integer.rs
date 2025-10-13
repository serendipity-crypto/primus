use primus_gcd::Xgcd;

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
    type SignedInteger: Integer;

    /// Returns `true` if and only if `self == 2^k` for some `k`.
    #[must_use]
    #[inline(always)]
    fn is_power_of_two(self) -> bool {
        self.count_ones() == 1
    }

    fn cast_from_signed(value: Self::SignedInteger) -> Self;

    /// Wrapping (modular) addition with a signed integer. Computes `self + rhs`, wrapping around at the boundary of the type.
    fn wrapping_add_signed(self, rhs: Self::SignedInteger) -> Self;
}

macro_rules! impl_unsigned_integer {
    ($t:ty, $i:ty) => {
        impl UnsignedInteger for $t {
            type SignedInteger = $i;

            #[inline]
            fn cast_from_signed(value: Self::SignedInteger) -> Self {
                value as $t
            }

            #[inline(always)]
            fn wrapping_add_signed(self, rhs: Self::SignedInteger) -> Self {
                <$t>::wrapping_add_signed(self, rhs)
            }
        }
    };
}

impl_unsigned_integer! {u8, i8}
impl_unsigned_integer! {u16, i16}
impl_unsigned_integer! {u32, i32}
impl_unsigned_integer! {u64, i64}
impl_unsigned_integer! {u128, i128}
