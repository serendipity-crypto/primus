use primus_integer::{CarryingMul, UnsignedInteger};

/// Computes `round(lhs * rhs / divisor)`.
#[inline]
pub(super) fn div_round<T>(lhs: T, rhs: T, divisor: T) -> T
where
    T: UnsignedInteger + CarryingMul,
{
    let (lo, hi) = lhs.carrying_mul(rhs, divisor >> 1u32);
    T::div_wide_fast(lo, hi, divisor)
}

/// Computes `round(lhs * rhs / divisor)` when the product fits in one limb.
#[inline]
pub(super) fn div_round_narrow<T>(lhs: T, rhs: T, divisor: T) -> T
where
    T: UnsignedInteger,
{
    let product = lhs * rhs;
    let (mut quotient, rem) = product.div_rem(divisor);
    if rem >= centered_half(divisor) {
        quotient += T::ONE;
    }
    quotient
}

#[inline]
pub(super) fn try_from_decoded<M, T>(decoded: T) -> M
where
    M: TryFrom<T>,
{
    M::try_from(decoded)
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}

#[inline]
pub(super) fn centered_half<T: UnsignedInteger>(t: T) -> T {
    (t >> 1u32) + (t & T::ONE)
}

#[inline]
pub(super) fn checked_message<M, T>(message: M) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}

#[inline]
pub(super) fn lift_centered<T>(message: T, t: T) -> (T, bool)
where
    T: UnsignedInteger,
{
    let half = centered_half(t);

    if message < half {
        (message, false)
    } else {
        (t - message, true)
    }
}

#[inline]
pub(super) fn lift_centered_from_raw<T: UnsignedInteger>(message: T, t: T, half: T) -> (T, bool) {
    if message < half {
        (message, false)
    } else {
        (t - message, true)
    }
}
