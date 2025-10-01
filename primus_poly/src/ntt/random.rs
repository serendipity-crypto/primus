use primus_integer::UnsignedInteger;
use primus_reduce::Modulus;
use rand::{CryptoRng, Rng, distr::Distribution};

use super::NttPolynomial;

impl<T: UnsignedInteger> NttPolynomial<T> {
    /// Generate a random [`NttPolynomial<T>`].
    #[inline]
    pub fn random<M, R>(n: usize, modulus: M, rng: &mut R) -> Self
    where
        M: Modulus<ValueT = T>,
        R: Rng + CryptoRng,
    {
        Self {
            values: modulus
                .uniform_distribution()
                .sample_iter(rng)
                .take(n)
                .collect(),
        }
    }

    /// Generate a random [`NttPolynomial<T>`]  with a specified distribution `distribution`.
    #[inline]
    pub fn random_with_distribution<R, D>(n: usize, rng: &mut R, distribution: D) -> Self
    where
        R: Rng + CryptoRng,
        D: Distribution<T>,
    {
        Self::new(distribution.sample_iter(rng).take(n).collect())
    }
}
