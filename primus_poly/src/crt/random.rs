use primus_distr::SignedDiscreteGaussian;
use primus_integer::{DataMut, RawData, UnsignedInteger};
use rand::distr::Uniform;

use super::CrtPolynomial;

impl<T: UnsignedInteger> CrtPolynomial<Vec<T>> {
    /// Generate a random binary [`CrtPolynomial<Vec<T>, T>`].
    #[inline]
    pub fn random_binary<R>(poly_length: usize, moduli_count: usize, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        Self(primus_distr::sample_crt_binary_values(
            poly_length,
            moduli_count,
            rng,
        ))
    }

    /// Generate a random ternary [`CrtPolynomial<Vec<T>, T>`].
    #[inline]
    pub fn random_ternary<R>(length: usize, moduli_minus_one: &[T], rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        Self(primus_distr::sample_crt_ternary_values(
            length,
            moduli_minus_one,
            rng,
        ))
    }

    /// Generate a random uniform [`CrtPolynomial<Vec<T>, T>`].
    #[inline]
    pub fn random_uniform<R>(length: usize, uniform_distrs: &[Uniform<T>], rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        Self(primus_distr::sample_crt_uniform_values(
            length,
            uniform_distrs,
            rng,
        ))
    }

    /// Generate a random gaussian [`CrtPolynomial<Vec<T>, T>`].
    #[inline]
    pub fn random_gaussian<R>(
        length: usize,
        moduli_value: &[T],
        gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        Self(primus_distr::sample_crt_gaussian_values(
            length,
            moduli_value,
            gaussian,
            rng,
        ))
    }
}

impl<S, T> CrtPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    #[inline]
    pub fn random_binary_assign<R>(&mut self, length: usize, rng: &mut R)
    where
        R: rand::Rng + rand::CryptoRng,
    {
        primus_distr::sample_crt_binary_values_inplace(self.as_mut(), length, rng)
    }

    #[inline]
    pub fn random_ternary_assign<R>(&mut self, length: usize, moduli_minus_one: &[T], rng: &mut R)
    where
        R: rand::Rng + rand::CryptoRng,
    {
        primus_distr::sample_crt_ternary_values_inplace(
            self.as_mut(),
            length,
            moduli_minus_one,
            rng,
        )
    }

    #[inline]
    pub fn random_uniform_assign<R>(
        &mut self,
        length: usize,
        uniform_distrs: &[Uniform<T>],
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
    {
        primus_distr::sample_crt_uniform_values_inplace(self.as_mut(), length, uniform_distrs, rng)
    }

    #[inline]
    pub fn random_gaussian_assign<R>(
        &mut self,
        length: usize,
        moduli_value: &[T],
        gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
    {
        primus_distr::sample_crt_gaussian_values_inplace(
            self.as_mut(),
            length,
            moduli_value,
            gaussian,
            rng,
        )
    }
}
