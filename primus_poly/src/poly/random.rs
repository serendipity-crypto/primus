use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_reduce::Modulus;
use rand::{CryptoRng, Rng, distr::Distribution};

use super::Polynomial;

impl<T: UnsignedInteger> Polynomial<T> {
    /// Generate a random [`Polynomial<T>`].
    #[inline]
    pub fn random<M, R>(n: usize, modulus: M, rng: &mut R) -> Self
    where
        M: Modulus<ValueT = T>,
        R: Rng + CryptoRng,
    {
        Self {
            poly: modulus
                .uniform_distribution()
                .sample_iter(rng)
                .take(n)
                .collect(),
        }
    }

    /// Generate a random [`Polynomial<T>`] with a specified `distribution`.
    #[inline]
    pub fn random_with_distribution<R, D>(n: usize, distribution: &D, rng: &mut R) -> Self
    where
        D: Distribution<T>,
        R: Rng + CryptoRng,
    {
        Self::new(distribution.sample_iter(rng).take(n).collect())
    }

    /// Generate a random binary [`Polynomial<T>`].
    #[inline]
    pub fn random_binary<R>(poly_length: usize, rng: &mut R) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self::new(primus_distr::sample_binary_values(poly_length, rng))
    }

    /// Generate a random ternary [`Polynomial<T>`].
    #[inline]
    pub fn random_ternary<R>(minus_one: T, poly_length: usize, rng: &mut R) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self::new(primus_distr::sample_ternary_values(
            minus_one,
            poly_length,
            rng,
        ))
    }

    /// Generate a random [`Polynomial<T>`] with discrete gaussian distribution.
    #[inline]
    pub fn random_gaussian<R>(
        poly_length: usize,
        gaussian: &DiscreteGaussian<T>,
        rng: &mut R,
    ) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self::new(gaussian.sample_iter(rng).take(poly_length).collect())
    }
}
