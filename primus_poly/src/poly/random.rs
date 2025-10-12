use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_reduce::{Modulus, ops::ReduceAddAssign};
use rand::{CryptoRng, Rng, distr::Distribution};

use crate::{ArrayBase, DataMut, DataOwned, RawData, poly::PolynomialOwned};

use super::Polynomial;

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Generate a random [`Polynomial<S, T>`].
    #[inline]
    pub fn random<M, R>(n: usize, modulus: M, rng: &mut R) -> Self
    where
        M: Modulus<ValueT = T>,
        R: Rng + CryptoRng,
    {
        Self(
            modulus
                .uniform_distribution()
                .sample_iter(rng)
                .take(n)
                .collect(),
        )
    }

    /// Generate a random [`Polynomial<S>`] with a specified `distribution`.
    #[inline]
    pub fn random_with_distribution<R, D>(n: usize, distribution: &D, rng: &mut R) -> Self
    where
        D: Distribution<T>,
        R: Rng + CryptoRng,
    {
        Self(distribution.sample_iter(rng).take(n).collect())
    }

    /// Generate a random [`Polynomial<S>`] with discrete gaussian distribution.
    #[inline]
    pub fn random_gaussian<R>(
        poly_length: usize,
        gaussian: &DiscreteGaussian<T>,
        rng: &mut R,
    ) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self(gaussian.sample_iter(rng).take(poly_length).collect())
    }
}

impl<T: UnsignedInteger> PolynomialOwned<T> {
    /// Generate a random binary [`Polynomial<S>`].
    #[inline]
    pub fn random_binary<R>(poly_length: usize, rng: &mut R) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self(ArrayBase(primus_distr::sample_binary_values(
            poly_length,
            rng,
        )))
    }

    /// Generate a random ternary [`Polynomial<S>`].
    #[inline]
    pub fn random_ternary<R>(minus_one: T, poly_length: usize, rng: &mut R) -> Self
    where
        R: Rng + CryptoRng,
    {
        Self(ArrayBase(primus_distr::sample_ternary_values(
            minus_one,
            poly_length,
            rng,
        )))
    }
}

impl<S, T> Polynomial<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Generate a random [`Polynomial<S, T>`].
    #[inline]
    pub fn random_assign<M, R>(&mut self, modulus: M, rng: &mut R)
    where
        M: Modulus<ValueT = T>,
        R: Rng + CryptoRng,
    {
        self.0
            .iter_mut()
            .zip(modulus.uniform_distribution().sample_iter(rng))
            .for_each(|(a, b)| *a = b);
    }

    /// Generate a random [`Polynomial<S>`] with a specified `distribution`.
    #[inline]
    pub fn random_with_distribution_assign<R, D>(&mut self, distribution: &D, rng: &mut R)
    where
        D: Distribution<T>,
        R: Rng + CryptoRng,
    {
        self.0
            .iter_mut()
            .zip(distribution.sample_iter(rng))
            .for_each(|(a, b)| *a = b);
    }

    /// Generate a random binary [`Polynomial<S>`].
    #[inline]
    pub fn random_binary_assign<R>(&mut self, rng: &mut R)
    where
        R: Rng + CryptoRng,
    {
        primus_distr::sample_binary_values_inplace(self.as_mut(), rng)
    }

    /// Generate a random ternary [`Polynomial<S>`].
    #[inline]
    pub fn random_ternary_assign<R>(&mut self, minus_one: T, rng: &mut R)
    where
        R: Rng + CryptoRng,
    {
        primus_distr::sample_ternary_values_inplace(self.as_mut(), minus_one, rng)
    }

    /// Generate a random [`Polynomial<S>`] with discrete gaussian distribution..
    #[inline]
    pub fn random_gaussian_assign<R>(&mut self, gaussian: &DiscreteGaussian<T>, rng: &mut R)
    where
        R: Rng + CryptoRng,
    {
        self.0
            .iter_mut()
            .zip(gaussian.sample_iter(rng))
            .for_each(|(a, b)| *a = b);
    }

    /// Generate a random [`Polynomial<S>`] with discrete gaussian distribution..
    #[inline]
    pub fn add_random_gaussian_assign<R, M>(
        &mut self,
        gaussian: &DiscreteGaussian<T>,
        modulus: M,
        rng: &mut R,
    ) where
        R: Rng + CryptoRng,
        M: Copy + ReduceAddAssign<T>,
    {
        self.0
            .iter_mut()
            .zip(gaussian.sample_iter(rng))
            .for_each(|(a, b)| modulus.reduce_add_assign(a, b));
    }
}
