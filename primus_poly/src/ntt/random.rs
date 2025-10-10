use primus_integer::UnsignedInteger;
use primus_reduce::Modulus;
use rand::{CryptoRng, Rng, distr::Distribution};

use crate::{DataOwned, RawData};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Generate a random [`NttPolynomial<T>`].
    #[inline]
    pub fn random<M, R>(poly_length: usize, modulus: M, rng: &mut R) -> Self
    where
        M: Modulus<ValueT = T>,
        R: Rng + CryptoRng,
    {
        Self(
            modulus
                .uniform_distribution()
                .sample_iter(rng)
                .take(poly_length)
                .collect(),
        )
    }

    /// Generate a random [`NttPolynomial<T>`]  with a specified distribution `distribution`.
    #[inline]
    pub fn random_with_distribution<R, D>(poly_length: usize, rng: &mut R, distribution: D) -> Self
    where
        R: Rng + CryptoRng,
        D: Distribution<T>,
    {
        Self(distribution.sample_iter(rng).take(poly_length).collect())
    }
}
