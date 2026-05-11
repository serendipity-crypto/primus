use core::fmt::Debug;

use primus_gcd::Xgcd;

use crate::{
    BorrowingSub, CarryingAdd, CarryingMul, DivRem, DivRemScalar, DivWide, Integer, WideningMul,
};

/// An abstraction over unsigned integer types.
///
/// `UnsignedInteger` extends [`Integer`] with operations that only make sense
/// for unsigned values: carrying add/sub, widening/carrying multiplication,
/// division with remainder, fast wide division, the extended GCD
/// ([`Xgcd`](primus_gcd::Xgcd)), and conversions from signed integers.
///
/// It is implemented for all standard Rust unsigned integer types (`u8`–`u128`,
/// `usize`) and serves as the principal value-type bound throughout the
/// crate hierarchy.
///
/// # Associated type
///
/// [`SignedInteger`](UnsignedInteger::SignedInteger) is the matching signed
/// type (e.g. `i64` for `u64`) and is used internally by algorithms such as
/// [`Xgcd::xgcd`](primus_gcd::Xgcd::xgcd) that need signed intermediate
/// cofactors.
pub trait UnsignedInteger:
    Integer
    + num_traits::Unsigned
    + CarryingAdd<CarryT = bool>
    + BorrowingSub<BorrowT = bool>
    + WideningMul
    + CarryingMul
    + DivRem
    + DivWide
    + DivRemScalar
    + Xgcd
    + TryFrom<usize, Error: Debug>
    + TryInto<usize, Error: Debug>
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

#[cfg(target_pointer_width = "64")]
impl UnsignedInteger for usize {
    type SignedInteger = i64;
    #[inline]
    fn cast_from_signed(value: Self::SignedInteger) -> Self {
        value as usize
    }
    #[inline(always)]
    fn wrapping_add_signed(self, rhs: Self::SignedInteger) -> Self {
        <usize>::wrapping_add_signed(self, rhs as isize)
    }
}

#[cfg(target_pointer_width = "32")]
impl UnsignedInteger for usize {
    type SignedInteger = i32;
    #[inline]
    fn cast_from_signed(value: Self::SignedInteger) -> Self {
        value as usize
    }
    #[inline(always)]
    fn wrapping_add_signed(self, rhs: Self::SignedInteger) -> Self {
        <usize>::wrapping_add_signed(self, rhs as isize)
    }
}
