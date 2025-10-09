use primus_integer::{AsInto, UnsignedInteger};

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
        None => encode_native(message, t),
        Some(q) if q.is_power_of_two() => encode_pow_of_2(message, t, q),
        Some(q) => encode_normal(message, t, q),
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
    debug_assert!(q.is_power_of_two() && t.is_power_of_two());
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

pub fn encode_normal<M, T>(message: M, t: T, q: T) -> T
where
    T: UnsignedInteger,
    M: TryInto<T>,
{
    let message: T = message
        .try_into()
        .map_err(|_| "out of range integral type conversion attempted")
        .unwrap();

    let q: f64 = q.as_into();
    let t: f64 = t.as_into();
    let m: f64 = message.as_into();

    (q / t * m).round().as_into()
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
        None => decode_native(cipher, t),
        Some(q) if q.is_power_of_two() => decode_pow_of_2(cipher, t, q),
        Some(q) => decode_normal(cipher, t, q),
    }
}

pub fn decode_normal<M, T>(cipher: T, t: T, q: T) -> M
where
    M: TryFrom<T>,
    T: UnsignedInteger,
{
    debug_assert!(t.is_power_of_two());
    let q_f: f64 = q.as_into();
    let t_f: f64 = t.as_into();
    let c: f64 = cipher.as_into();
    let temp: T = (c / (q_f / t_f)).round().as_into();
    let temp = if temp >= t { temp - t } else { temp };

    M::try_from(temp)
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
    debug_assert!(q.is_power_of_two() && t.is_power_of_two());
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
