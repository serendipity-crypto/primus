use primus_integer::Integer;
use rand::distr::Distribution;

use crate::DistrErr;

mod cdt;
#[cfg(target_os = "linux")]
mod unix_cdt;
mod ziggurat;

/// The gaussian distribution `N(mean, std_dev**2)`.
#[derive(Clone)]
pub enum SignedDiscreteGaussian<T: Integer> {
    /// CDTSampler
    Cdt(cdt::CDTSampler<T>),
    /// DiscreteZiggurat
    Ziggurat(ziggurat::DiscreteZiggurat<T>),
}

impl<T: Integer> SignedDiscreteGaussian<T> {
    /// Construct, from mean and standard deviation
    ///
    /// Parameters:
    ///
    /// -   mean (`μ`, unrestricted)
    /// -   standard deviation (`σ`, must be finite)
    #[inline]
    pub fn new(std_dev: f64) -> Result<SignedDiscreteGaussian<T>, DistrErr<T>> {
        if std_dev < 16.0 {
            Ok(SignedDiscreteGaussian::Cdt(cdt::CDTSampler::new(
                std_dev, 12.0,
            )))
        } else {
            Ok(SignedDiscreteGaussian::Ziggurat(
                ziggurat::DiscreteZiggurat::new(std_dev, 12.0),
            ))
        }
    }

    /// Construct, from mean and standard deviation
    ///
    /// Parameters:
    ///
    /// -   mean (`μ`, unrestricted)
    /// -   standard deviation (`σ`, must be finite)
    #[inline]
    pub fn new_with_max_limit(
        std_dev: f64,
        max_std_dev: f64,
    ) -> Result<SignedDiscreteGaussian<T>, DistrErr<T>> {
        if max_std_dev <= std_dev || std_dev < 0.7 {
            return Err(DistrErr::DiscreteGaussianErr {
                std_dev,
                modulus_minus_one: T::ZERO,
            });
        }
        unimplemented!()
    }

    /// Returns the standard deviation of this [`SignedDiscreteGaussian<T>`].
    pub fn standard_deviation(&self) -> f64 {
        match self {
            SignedDiscreteGaussian::Cdt(sampler) => sampler.std_dev(),
            SignedDiscreteGaussian::Ziggurat(sampler) => sampler.std_dev(),
        }
    }
}

impl<T: Integer> Distribution<T> for SignedDiscreteGaussian<T> {
    #[inline]
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        match self {
            SignedDiscreteGaussian::Cdt(sampler) => sampler.sample(rng),
            SignedDiscreteGaussian::Ziggurat(sampler) => sampler.sample(rng),
        }
    }
}
