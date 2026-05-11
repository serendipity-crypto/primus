//! SIMD abstractions for unsigned integer types.
//!
//! This module provides traits and blanket implementations that extend the
//! scalar [`UnsignedInteger`](crate::UnsignedInteger) operations to SIMD
//! vectors when the `nightly` and `simd` features are enabled.
//!
//! [`SimdUnsignedInteger`] marks unsigned integer types that can serve as
//! SIMD lane elements. [`SimdArray`] extends [`Simd`] vectors with the
//! arithmetic and comparison capabilities required by higher-level crates.
//! [`SimdMaskArray`] provides the corresponding mask operations.

use core::{
    fmt::Debug,
    iter::{Product, Sum},
    ops::*,
    simd::{
        LaneCount, Mask, MaskElement, Simd, SimdCast, SimdElement, SupportedLaneCount,
        cmp::{SimdOrd, SimdPartialEq, SimdPartialOrd},
        num::SimdUint,
    },
};

use crate::{BorrowingSub, CarryingAdd, CarryingMul, WideningMul};

use super::UnsignedInteger;

pub trait SimdUnsignedInteger: UnsignedInteger + SimdElement + SimdCast {}

macro_rules! impl_simd_unsigned_integer {
    ($($t:ty)*) => ($(
        impl SimdUnsignedInteger for $t {}
    )*)
}

impl_simd_unsigned_integer! {u8 u16 u32 u64 usize}

pub trait SimdArray<T: SimdUnsignedInteger, const N: usize>
where
    LaneCount<N>: SupportedLaneCount,
    Self: Send + Sync + Clone + Copy + Default,
    Self: PartialEq + PartialOrd + Eq + Ord,
    Self: Debug,
    Self: From<Simd<T, N>> + Into<Simd<T, N>>,
    Self: From<[T; N]> + Into<[T; N]>,
    Self: AsRef<[T; N]> + AsMut<[T; N]>,
    Self:
        SimdPartialEq<Mask: SimdMaskArray<<T as SimdElement>::Mask, N>> + SimdPartialOrd + SimdOrd,
    Self: Product<Self> + Sum<Self>,
    for<'a> Self: Product<&'a Self> + Sum<&'a Self>,
    Self: Index<usize, Output = T> + IndexMut<usize, Output = T>,
    Self: Add<Output = Self> + AddAssign,
    Self: Sub<Output = Self> + SubAssign,
    Self: Mul<Output = Self> + MulAssign,
    Self: Div<Output = Self> + DivAssign,
    Self: Rem<Output = Self> + RemAssign,
    for<'a> Self: Add<&'a Self, Output = Self> + AddAssign<&'a Self>,
    for<'a> Self: Sub<&'a Self, Output = Self> + SubAssign<&'a Self>,
    for<'a> Self: Mul<&'a Self, Output = Self> + MulAssign<&'a Self>,
    for<'a> Self: Div<&'a Self, Output = Self> + DivAssign<&'a Self>,
    for<'a> Self: Rem<&'a Self, Output = Self> + RemAssign<&'a Self>,
    Self: SimdUint,
    Self: CarryingAdd<CarryT = Self::Mask>
        + BorrowingSub<BorrowT = Self::Mask>
        + WideningMul
        + CarryingMul,
{
    /// Returns a tuple of the sum along with a boolean indicating whether an arithmetic overflow would occur.
    /// If an overflow would have occurred then the wrapped value is returned.
    fn overflowing_add(self, rhs: Self) -> (Self, Self::Mask) {
        let a = self + rhs;
        (a, a.simd_lt(self))
    }
}

macro_rules! impl_simd_array {
    ($($t:ty)*) => ($(
        impl<const N: usize> SimdArray<$t, N> for Simd<$t, N> where LaneCount<N>: SupportedLaneCount {}
    )*)
}

impl_simd_array! {u8 u16 u32 u64 usize}

#[allow(clippy::len_without_is_empty)]
pub trait SimdMaskArray<T: MaskElement, const N: usize>
where
    LaneCount<N>: SupportedLaneCount,
    Self: Send + Sync + Clone + Copy + Default,
    Self: PartialEq + PartialOrd,
    Self: Debug,
    Self: From<[bool; N]> + Into<[bool; N]>,
    Self: SimdPartialEq<Mask = Self> + SimdPartialOrd + SimdOrd,
    Self: BitAnd<Output = Self> + BitAndAssign + BitAnd<bool, Output = Self> + BitAndAssign<bool>,
    Self: BitOr<Output = Self> + BitOrAssign + BitOr<bool, Output = Self> + BitOrAssign<bool>,
    Self: BitXor<Output = Self> + BitXorAssign + BitXor<bool, Output = Self> + BitXorAssign<bool>,
    Self: Not<Output = Self>,
{
    /// Get the number of lanes in this vector.
    #[must_use]
    #[inline]
    fn len(&self) -> usize {
        N
    }

    /// Choose elements from two vectors.
    ///
    /// For each element in the mask, choose the corresponding element from `true_values` if
    /// that element mask is true, and `false_values` if that element mask is false.
    ///
    /// # Examples
    /// ```ignore
    /// # #![feature(portable_simd)]
    /// # use core::simd::{Simd, Mask};
    /// let a = Simd::from_array([0, 1, 2, 3]);
    /// let b = Simd::from_array([4, 5, 6, 7]);
    /// let mask = Mask::from_array([true, false, false, true]);
    /// let c = mask.select(a, b);
    /// assert_eq!(c.to_array(), [0, 5, 6, 3]);
    /// ```
    #[must_use = "method returns a new vector and does not mutate the original inputs"]
    fn select<U>(self, true_values: Simd<U, N>, false_values: Simd<U, N>) -> Simd<U, N>
    where
        U: SimdElement<Mask = T>;

    /// Constructs a mask by setting all elements to the given value.
    fn splat(value: bool) -> Self;

    /// Converts an array of bools to a SIMD mask.
    fn from_array(array: [bool; N]) -> Self;

    /// Converts a SIMD mask to an array of bools.
    fn to_array(self) -> [bool; N];

    /// Converts a vector of integers to a mask, where 0 represents `false` and -1
    /// represents `true`.
    ///
    /// # Panics
    /// Panics if any element is not 0 or -1.
    #[must_use = "method returns a new mask and does not mutate the original value"]
    #[track_caller]
    fn from_int(value: Simd<T, N>) -> Self;

    /// Converts the mask to a vector of integers, where 0 represents `false` and -1
    /// represents `true`.
    #[must_use = "method returns a new vector and does not mutate the original value"]
    fn to_int(self) -> Simd<T, N>;

    /// Returns true if any element is set, or false otherwise.
    #[must_use = "method returns a new bool and does not mutate the original value"]
    fn any(self) -> bool;

    /// Returns true if all elements are set, or false otherwise.
    #[must_use = "method returns a new bool and does not mutate the original value"]
    fn all(self) -> bool;
}

macro_rules! impl_mask_array {
    ($($t:ty)*) => ($(
        impl<const N: usize> SimdMaskArray<$t, N> for Mask<$t, N>
        where
            LaneCount<N>: SupportedLaneCount
        {
            #[inline]
            fn select<U>(self, true_values: Simd<U, N>, false_values: Simd<U, N>) -> Simd<U, N>
            where
                U: SimdElement<Mask = $t>
            {
                self.select(true_values, false_values)
            }

            #[inline]
            fn splat(value: bool) -> Self {
                Self::splat(value)
            }

            #[inline]
            fn from_array(array: [bool; N]) -> Self {
                Self::from_array(array)
            }

            #[inline]
            fn to_array(self) -> [bool; N] {
                self.to_array()
            }

            #[inline]
            fn from_int(value: Simd<$t, N>) -> Self {
                Self::from_int(value)
            }

            #[inline]
            fn to_int(self) -> Simd<$t, N> {
                self.to_int()
            }

            #[inline]
            fn any(self) -> bool {
                self.any()
            }

            #[inline]
            fn all(self) -> bool {
                self.all()
            }
        }
    )*)
}

impl_mask_array! {i8 i16 i32 i64 isize}
