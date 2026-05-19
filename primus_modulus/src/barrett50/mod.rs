//! Barrett reduction specialized for 50-bit moduli with an optional
//! AVX-512-IFMA fast path.
//!
//! Construct with [`Barrett50Modulus::new`] requiring `m ∈ [2^48, 2^50)`.
//! Slice operations follow the trait surface of [`super::BarrettModulus`]; on
//! hosts compiled with `target_feature = "avx512ifma"` (and `avx512f` /
//! `avx512dq`) the mul-family slice kernels go through the IFMA path in
//! [`simd_ifma`]. Other configurations transparently delegate to the inner
//! [`super::BarrettModulus<u64>`].

use core::fmt::Display;

use primus_integer::DivRemScalar;

use crate::BarrettModulus;

mod ops;
mod slice;

#[cfg(all(
    feature = "nightly",
    feature = "simd",
    target_feature = "avx512f",
    target_feature = "avx512dq",
    target_feature = "avx512ifma",
))]
mod simd_ifma;

/// Inclusive lower bound on the modulus: keeps `µ_hi ≤ 16` so
/// `x_hi · µ_hi` fits in 52 bits.
pub const MIN_VALUE: u64 = 1u64 << 48;
/// Exclusive upper bound on the modulus: keeps `a · b < 2^100` so the
/// IFMA `madd52lo`/`madd52hi` pair captures the full product.
pub const MAX_VALUE: u64 = 1u64 << 50;

/// Barrett reduction wrapper for moduli `m ∈ [2^48, 2^50)`.
///
/// Stores two precomputed forms of `µ = ⌊2^104 / m⌋`:
///
/// - `mu_lo52` / `mu_hi`: 52-bit / ≤4-bit pieces consumed by the IFMA
///   kernels in [`simd_ifma`].
/// - `inner`: a full [`BarrettModulus<u64>`] reused for scalar ops and as
///   the non-IFMA SIMD fallback.
#[derive(Debug, Clone, Copy)]
pub struct Barrett50Modulus {
    pub(super) value: u64,
    pub(super) inner: BarrettModulus<u64>,
    /// Low 52 bits of `⌊2^104 / m⌋`.
    #[allow(dead_code, reason = "consumed by IFMA path; unused in non-IFMA builds")]
    pub(super) mu_lo52: u64,
    /// `⌊2^104 / m⌋ >> 52`. For `m ∈ [2^48, 2^50)` this is in `[4, 16]`.
    #[allow(dead_code, reason = "consumed by IFMA path; unused in non-IFMA builds")]
    pub(super) mu_hi: u64,
    /// `(1 << 52) − m`, used by `madd52lo` to implement subtraction.
    #[allow(dead_code, reason = "consumed by IFMA path; unused in non-IFMA builds")]
    pub(super) neg_m_mod_pow2_52: u64,
}

impl Barrett50Modulus {
    /// Construct a [`Barrett50Modulus`] for `value ∈ [2^48, 2^50)`.
    ///
    /// # Panics
    ///
    /// Panics if `value < 2^48` or `value >= 2^50`.
    #[inline]
    pub fn new(value: u64) -> Self {
        assert!(
            (MIN_VALUE..MAX_VALUE).contains(&value),
            "Barrett50Modulus requires value in [2^48, 2^50); got {value}"
        );
        Self::new_unchecked(value)
    }

    /// Fallible constructor returning `None` if the value is out of range.
    #[inline]
    pub fn try_new(value: u64) -> Option<Self> {
        if !(MIN_VALUE..MAX_VALUE).contains(&value) {
            return None;
        }
        Some(Self::new_unchecked(value))
    }

    /// Like [`new`](Self::new) but without the range check.
    ///
    /// # Correctness
    ///
    /// `value` must satisfy `2^48 ≤ value < 2^50`. Violating this breaks
    /// the IFMA path: products may exceed 2^100, and `µ_hi` may overflow
    /// `u64` lane multiplications in the slice kernels.
    #[inline]
    pub fn new_unchecked(value: u64) -> Self {
        // µ_104 = ⌊2^104 / value⌋. Dividend laid out as little-endian
        // 64-bit limbs: 2^104 = 0 + (1 << 40) · 2^64 + 0 · 2^128.
        let mut quotient = [0u64; 3];
        let _rem =
            <u64 as DivRemScalar>::div_rem_scalar(&[0u64, 1u64 << 40, 0u64], value, &mut quotient);
        // For value ∈ [2^48, 2^50): µ_104 ∈ (2^54, 2^56], so quotient[1..]
        // is all zero — the value fits in quotient[0].
        let mu_104 = quotient[0];
        let mu_lo52 = mu_104 & ((1u64 << 52) - 1);
        let mu_hi = mu_104 >> 52;
        let neg_m_mod_pow2_52 = (1u64 << 52).wrapping_sub(value);

        Self {
            value,
            inner: BarrettModulus::<u64>::new(value),
            mu_lo52,
            mu_hi,
            neg_m_mod_pow2_52,
        }
    }

    /// Returns the modulus.
    #[inline]
    pub const fn value(&self) -> u64 {
        self.value
    }

    /// Returns the wrapped [`BarrettModulus<u64>`].
    #[inline]
    pub const fn inner(&self) -> BarrettModulus<u64> {
        self.inner
    }
}

impl From<Barrett50Modulus> for BarrettModulus<u64> {
    #[inline]
    fn from(m: Barrett50Modulus) -> Self {
        m.inner
    }
}

impl Display for Barrett50Modulus {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl primus_reduce::Modulus for Barrett50Modulus {
    type ValueT = u64;

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
        self.value - 1
    }
}

#[cfg(test)]
mod tests {
    use primus_reduce::FieldContext;

    use super::*;

    fn field_trait<T: primus_integer::UnsignedInteger, M: FieldContext<T>>(_modulus: M) {}

    #[test]
    fn implements_field_context() {
        field_trait(Barrett50Modulus::new(1u64 << 48));
        field_trait(Barrett50Modulus::new((1u64 << 50) - 27));
    }

    #[test]
    fn precompute_matches_2_to_104_div_value() {
        for &value in &[
            1u64 << 48,
            (1u64 << 48) + 1,
            (1u64 << 49) + 17,
            (1u64 << 50) - 27,
            (1u64 << 50) - 1,
        ] {
            let m = Barrett50Modulus::new(value);
            let mu = (m.mu_hi as u128) << 52 | m.mu_lo52 as u128;
            let expected = (1u128 << 104) / value as u128;
            assert_eq!(mu, expected, "value={value}");
            assert_eq!(m.neg_m_mod_pow2_52, (1u64 << 52) - value);
            assert!(m.mu_hi <= 16);
        }
    }

    #[test]
    #[should_panic]
    fn rejects_below_2_to_48() {
        let _ = Barrett50Modulus::new((1u64 << 48) - 1);
    }

    #[test]
    #[should_panic]
    fn rejects_2_to_50() {
        let _ = Barrett50Modulus::new(1u64 << 50);
    }

    #[test]
    fn try_new_returns_none_out_of_range() {
        assert!(Barrett50Modulus::try_new((1u64 << 48) - 1).is_none());
        assert!(Barrett50Modulus::try_new(1u64 << 50).is_none());
        assert!(Barrett50Modulus::try_new(1u64 << 48).is_some());
        assert!(Barrett50Modulus::try_new((1u64 << 50) - 1).is_some());
    }
}
