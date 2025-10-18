use std::slice::IterMut;

use num_traits::ConstZero;
use primus_integer::{AsInto, Integer, UnsignedInteger};
use rand::{
    CryptoRng, Rng,
    distr::{Distribution, Uniform},
};

use crate::SignedDiscreteGaussian;

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

pub fn sample_crt_binary_values<T, R>(length: usize, moduli_count: usize, rng: &mut R) -> Vec<T>
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let mut result = vec![T::ZERO; length * moduli_count];

    let (v, w) = result.split_at_mut(length);

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

    w.chunks_exact_mut(length)
        .for_each(|s| s.copy_from_slice(v));

    result
}

pub fn sample_crt_binary_values_inplace<T, R>(result: &mut [T], chunk_size: usize, rng: &mut R)
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let (v, w) = result.split_at_mut(chunk_size);

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

    w.chunks_exact_mut(chunk_size)
        .for_each(|s| s.copy_from_slice(v));
}

/// Sample a ternary vector whose values are `T`.
pub fn sample_crt_ternary_values<T, R>(length: usize, moduli_minus_one: &[T], rng: &mut R) -> Vec<T>
where
    T: Integer,
    R: Rng + CryptoRng,
{
    let moduli_count = moduli_minus_one.len();
    let mut result = vec![T::ZERO; length * moduli_count];

    let mut iters: Vec<IterMut<'_, T>> = result
        .chunks_exact_mut(length)
        .map(|s| s.iter_mut())
        .collect();

    'outer: loop {
        let mut r = rng.next_u32();
        for _ in 0..16 {
            match r & 0b11 {
                0 | 1 => {
                    for iter in iters.iter_mut() {
                        if let Some(value) = iter.next() {
                            *value = T::ZERO;
                        } else {
                            break 'outer;
                        }
                    }
                }
                2 => {
                    for iter in iters.iter_mut() {
                        if let Some(value) = iter.next() {
                            *value = T::ONE;
                        } else {
                            break 'outer;
                        }
                    }
                }
                3 => {
                    for (i, iter) in iters.iter_mut().enumerate() {
                        if let Some(value) = iter.next() {
                            *value = moduli_minus_one[i];
                        } else {
                            break 'outer;
                        }
                    }
                }
                _ => unreachable!(),
            }
            r >>= 2;
        }
    }

    result
}

/// Sample a gaussian vector whose values are `T`.
pub fn sample_crt_gaussian_values<T, R>(
    length: usize,
    moduli: &[T],
    gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
    rng: &mut R,
) -> Vec<T>
where
    T: UnsignedInteger,
    R: Rng + CryptoRng,
{
    let bound: f64 = 24.0 * gaussian.standard_deviation();
    let bound: T = bound.as_into();
    for modulus in moduli {
        assert!(bound < *modulus);
    }

    let moduli_count = moduli.len();
    let mut result = vec![T::ZERO; length * moduli_count];

    let mut iters: Vec<IterMut<'_, T>> = result
        .chunks_exact_mut(length)
        .map(|s| s.iter_mut())
        .collect();

    'outer: loop {
        let r = gaussian.sample(rng);
        if r >= <<T as UnsignedInteger>::SignedInteger as ConstZero>::ZERO {
            let t: T = T::cast_from_signed(r);
            for iter in iters.iter_mut() {
                if let Some(value) = iter.next() {
                    *value = t;
                } else {
                    break 'outer;
                }
            }
        } else {
            for (iter, modulus) in iters.iter_mut().zip(moduli) {
                if let Some(value) = iter.next() {
                    *value = modulus.wrapping_add_signed(r);
                } else {
                    break 'outer;
                }
            }
        }
    }

    result
}

/// Sample a gaussian vector whose values are `T`.
pub fn sample_crt_gaussian_values_inplace<T, R>(
    result: &mut [T],
    chunk_size: usize,
    moduli: &[T],
    gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
    rng: &mut R,
) where
    T: UnsignedInteger,
    R: Rng + CryptoRng,
{
    let iters: Vec<IterMut<'_, T>> = result
        .chunks_exact_mut(chunk_size)
        .map(|s| s.iter_mut())
        .collect();
    sample_crt_gaussian_values_iter_mut(iters, moduli, gaussian, rng);
}

/// Sample a gaussian vector whose values are `T`.
pub fn sample_crt_gaussian_values_iter_mut<T, R>(
    mut iters: Vec<IterMut<'_, T>>,
    moduli: &[T],
    gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
    rng: &mut R,
) where
    T: UnsignedInteger,
    R: Rng + CryptoRng,
{
    loop {
        let r = gaussian.sample(rng);
        if r >= <<T as UnsignedInteger>::SignedInteger as ConstZero>::ZERO {
            let t: T = T::cast_from_signed(r);
            for iter in iters.iter_mut() {
                if let Some(value) = iter.next() {
                    *value = t;
                } else {
                    return;
                }
            }
        } else {
            for (iter, modulus) in iters.iter_mut().zip(moduli) {
                if let Some(value) = iter.next() {
                    *value = modulus.wrapping_add_signed(r);
                } else {
                    return;
                }
            }
        }
    }
}

/// Sample a uniform vector whose values are `T`.
pub fn sample_crt_uniform_values_iter_mut<T, R>(
    iters: Vec<IterMut<'_, T>>,
    uniform_distrs: &[Uniform<T>],
    rng: &mut R,
) where
    T: UnsignedInteger,
    R: Rng + CryptoRng,
{
    iters.into_iter().zip(uniform_distrs).for_each(|(s, u)| {
        s.zip(u.sample_iter(&mut *rng)).for_each(|(a, b)| {
            *a = b;
        });
    });
}
