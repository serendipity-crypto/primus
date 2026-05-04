use primus_integer::UnsignedInteger;

#[inline]
fn div_wide<T: UnsignedInteger>(lo: T, hi: T, divisor: T) -> T {
    let mut quotient = [T::ZERO; 2];
    T::div_rem_scalar(&[lo, hi], divisor, &mut quotient);
    quotient[0]
}

#[inline]
fn mul_div_round<T: UnsignedInteger>(lhs: T, rhs: T, divisor: T) -> T {
    let (lo, hi) = lhs.carrying_mul(rhs, divisor >> 1u32);
    div_wide(lo, hi, divisor)
}

/// Encodes a message.
///
/// # Parameters
///
/// - `t` is message space
/// - `q` is LWE modulus value.
#[inline]
pub fn encode<M, T>(message: M, t: T, q: Option<T>) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    match q {
        None if t.is_power_of_two() => encode_native(message, t),
        None => encode_native_scaled(message, t),
        Some(q) if q.is_power_of_two() && t.is_power_of_two() => encode_pow_of_2(message, t, q),
        Some(q) => encode_scaled(message, t, q),
    }
}

/// Encodes a message.
///
/// # Parameters
///
/// - `t` is message space
/// - `q` is LWE modulus value.
/// - This function needs `q` and `t` are power of 2.
///
/// # Panic
///
/// Panics if the message exceeds the message space.
#[inline]
pub fn encode_pow_of_2<M, T>(message: M, t: T, q: T) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    debug_assert!(q.is_power_of_two() && t.is_power_of_two() && t < q);
    // Shift the message to the most significant part of `T`.
    let message: T = message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap();
    assert!(
        message < t,
        "message {message} is bigger than the message space"
    );
    message << (q / t).trailing_zeros()
}

/// Encodes a message.
///
/// # Parameters
///
/// - `t` is message space
/// - This function needs `t` be power of 2.
///
/// # Panic
///
/// Panics if the message exceeds the message space.
#[inline]
pub fn encode_native<M, T>(message: M, t: T) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    debug_assert!(t.is_power_of_two());
    let message: T = message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap();
    assert!(
        message < t,
        "message {message} is bigger than the message space"
    );
    message << (T::BITS - t.trailing_zeros())
}

/// Encodes a message under the native modulus `2^T::BITS` for arbitrary `t`.
///
/// Computes `round(message * 2^T::BITS / t)` using integer arithmetic.
///
/// # Parameters
///
/// - `message` is the plaintext message.
/// - `t` is the message space.
///
/// # Panic
///
/// Panics if the message exceeds the message space.
#[inline]
pub fn encode_native_scaled<M, T>(message: M, t: T) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    let message: T = message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap();
    assert!(message < t);

    div_wide(t >> 1u32, message, t)
}

/// Encodes a message by scaling it from the message space into modulus `q`.
///
/// Computes `round(message * q / t)` using integer arithmetic.
///
/// # Parameters
///
/// - `message` is the plaintext message.
/// - `t` is the message space.
/// - `q` is the ciphertext modulus value.
///
/// # Panic
///
/// Panics if `t` or `q` is zero, or if the message exceeds the message space.
#[inline]
pub fn encode_scaled<M, T>(message: M, t: T, q: T) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    assert!(!t.is_zero());
    assert!(!q.is_zero());

    let message: T = message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap();
    assert!(message < t);

    mul_div_round(message, q, t)
}

/// Decodes an encode value.
///
/// # Parameters
///
/// - `t` is message space
/// - `q` is LWE modulus value.
#[inline]
pub fn decode<M, T>(cipher: T, t: T, q: Option<T>) -> M
where
    M: TryFrom<T>,
    T: UnsignedInteger,
{
    match q {
        None if t.is_power_of_two() => decode_native(cipher, t),
        None => decode_native_scaled(cipher, t),
        Some(q) if q.is_power_of_two() && t.is_power_of_two() => decode_pow_of_2(cipher, t, q),
        Some(q) => decode_scaled(cipher, t, q),
    }
}

/// Decodes a scaled plaintext value from modulus `q` into the message space.
///
/// Computes `round(cipher * t / q) mod t` using integer arithmetic.
///
/// # Parameters
///
/// - `cipher` is the encoded plaintext value under modulus `q`.
/// - `t` is the message space.
/// - `q` is the ciphertext modulus value.
///
/// # Panic
///
/// Panics if the decoded message cannot fit in `M`.
#[inline]
pub fn decode_scaled<M, T>(cipher: T, t: T, q: T) -> M
where
    T: UnsignedInteger,
    M: TryFrom<T>,
{
    let mut decoded = mul_div_round(cipher, t, q);

    if decoded >= t {
        decoded -= t;
    }

    M::try_from(decoded)
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}

/// Decodes an encode value.
///
/// # Parameters
///
/// - `t` is message space
/// - `q` is LWE modulus value.
/// - This function needs `q` and `t` are power of 2.
///
/// # Panic
///
/// Panics if the decoded message cannot fit in `M`.
#[inline]
pub fn decode_pow_of_2<M, T>(cipher: T, t: T, q: T) -> M
where
    M: TryFrom<T>,
    T: UnsignedInteger,
{
    debug_assert!(q.is_power_of_two() && t.is_power_of_two() && t < q);
    // Move the message to the least significant part of `C`.
    // Leave one more bit for round.
    let temp = cipher >> ((q / t).trailing_zeros() - 1);
    let decoded = ((temp + T::ONE) >> 1u32) & (t - T::ONE);

    M::try_from(decoded)
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}

/// Decodes an encode value.
///
/// # Parameters
///
/// - `t` is message space
/// - `q` is LWE modulus value.
/// - This function needs `t` be power of 2.
///
/// # Panic
///
/// Panics if the decoded message cannot fit in `M`.
#[inline]
pub fn decode_native<M, T>(cipher: T, t: T) -> M
where
    M: TryFrom<T>,
    T: UnsignedInteger,
{
    debug_assert!(t.is_power_of_two());
    // Move the message to the least significant part of `C`.
    // Leave one more bit for round.
    let temp = cipher >> (T::BITS - t.trailing_zeros() - 1);
    let decoded = ((temp + T::ONE) >> 1u32) & (t - T::ONE);

    M::try_from(decoded)
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}

/// Decodes a scaled plaintext value from the native modulus `2^T::BITS`.
///
/// Computes `round(cipher * t / 2^T::BITS) mod t` using integer arithmetic.
///
/// # Parameters
///
/// - `cipher` is the encoded plaintext value under the native modulus.
/// - `t` is the message space.
///
/// # Panic
///
/// Panics if the decoded message cannot fit in `M`.
#[inline]
pub fn decode_native_scaled<M, T>(cipher: T, t: T) -> M
where
    T: UnsignedInteger,
    M: TryFrom<T>,
{
    let mut decoded = cipher.carrying_mul_hw(t, T::ONE << (T::BITS - 1));

    if decoded >= t {
        decoded -= t;
    }

    M::try_from(decoded)
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap()
}
