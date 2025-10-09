use primus_integer::Integer;
use rand::{CryptoRng, Rng};

/// Sample a binary vector whose values are `T`.
pub fn sample_binary_values<T, R>(length: usize, rng: &mut R) -> Vec<T>
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let mut v = vec![T::ZERO; length];
    let mut iter = v.chunks_exact_mut(32);
    for chunk in &mut iter {
        let mut r = rng.next_u32();
        for elem in chunk.iter_mut() {
            *elem = T::as_from(r & 0b1);
            r >>= 1;
        }
    }
    let mut r = rng.next_u32();
    for elem in iter.into_remainder() {
        *elem = T::as_from(r & 0b1);
        r >>= 1;
    }
    v
}

/// Sample a binary vector whose values are `T`.
pub fn sample_binary_values_inplace<T, R>(result: &mut [T], rng: &mut R)
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let mut iter = result.chunks_exact_mut(32);
    for chunk in &mut iter {
        let mut r = rng.next_u32();
        for elem in chunk.iter_mut() {
            *elem = T::as_from(r & 0b1);
            r >>= 1;
        }
    }
    let mut r = rng.next_u32();
    for elem in iter.into_remainder() {
        *elem = T::as_from(r & 0b1);
        r >>= 1;
    }
}

/// Sample a ternary vector whose values are `T`.
pub fn sample_ternary_values<T, R>(minus_one: T, length: usize, rng: &mut R) -> Vec<T>
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let s = [T::ZERO, T::ZERO, T::ONE, minus_one];
    let mut v = vec![T::ZERO; length];
    let mut iter = v.chunks_exact_mut(16);
    for chunk in &mut iter {
        let mut r = rng.next_u32();
        for elem in chunk.iter_mut() {
            *elem = s[(r & 0b11) as usize];
            r >>= 2;
        }
    }
    let mut r = rng.next_u32();
    for elem in iter.into_remainder() {
        *elem = s[(r & 0b11) as usize];
        r >>= 2;
    }
    v
}

/// Sample a ternary vector whose values are `T`.
pub fn sample_ternary_values_inplace<T, R>(result: &mut [T], minus_one: T, rng: &mut R)
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let s = [T::ZERO, T::ZERO, T::ONE, minus_one];

    let mut iter = result.chunks_exact_mut(16);
    for chunk in &mut iter {
        let mut r = rng.next_u32();
        for elem in chunk.iter_mut() {
            *elem = s[(r & 0b11) as usize];
            r >>= 2;
        }
    }
    let mut r = rng.next_u32();
    for elem in iter.into_remainder() {
        *elem = s[(r & 0b11) as usize];
        r >>= 2;
    }
}
