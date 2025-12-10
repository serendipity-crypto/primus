use primus_integer::UnsignedInteger;
use rand::{Rng, distr::Distribution};

mod cdt;
#[cfg(all(target_os = "linux", feature = "high_precision"))]
mod unix_cdt;
mod ziggurat;

pub use cdt::CDTSampler;
#[cfg(all(target_os = "linux", feature = "high_precision"))]
pub use unix_cdt::UnixCDTSampler;
pub use ziggurat::DiscreteZiggurat;

use crate::DistrErr;

/// The gaussian distribution `N(mean, std_dev**2)`.
#[derive(Clone)]
pub enum DiscreteGaussian<T: UnsignedInteger> {
    /// CDTSampler
    Cdt(CDTSampler<T>),
    /// DiscreteZiggurat
    Ziggurat(DiscreteZiggurat<T>),
}

impl<T: UnsignedInteger> DiscreteGaussian<T> {
    /// Construct, from mean and standard deviation
    ///
    /// Parameters:
    ///
    /// -   mean (`μ`, unrestricted)
    /// -   standard deviation (`σ`, must be finite)
    #[inline]
    pub fn new(std_dev: f64, modulus_minus_one: T) -> Result<DiscreteGaussian<T>, DistrErr<T>> {
        if std_dev < 0.7 {
            Err(DistrErr::DiscreteGaussianErr {
                std_dev,
                modulus_minus_one,
            })
        } else if std_dev <= 16.0 {
            Ok(DiscreteGaussian::Cdt(CDTSampler::new(
                std_dev,
                12.0,
                modulus_minus_one,
            )))
        } else {
            Ok(DiscreteGaussian::Ziggurat(DiscreteZiggurat::new(
                std_dev,
                12.0,
                modulus_minus_one,
            )))
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
        modulus_minus_one: T,
    ) -> Result<DiscreteGaussian<T>, DistrErr<T>> {
        if max_std_dev <= std_dev || std_dev < 0.7 {
            return Err(DistrErr::DiscreteGaussianErr {
                std_dev,
                modulus_minus_one,
            });
        }
        unimplemented!()
    }

    /// Returns the standard deviation of this [`DiscreteGaussian<T>`].
    pub fn standard_deviation(&self) -> f64 {
        match self {
            DiscreteGaussian::Cdt(sampler) => sampler.std_dev(),
            DiscreteGaussian::Ziggurat(sampler) => sampler.std_dev(),
        }
    }
}

impl<T: UnsignedInteger> Distribution<T> for DiscreteGaussian<T> {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> T {
        match self {
            DiscreteGaussian::Cdt(sampler) => sampler.sample(rng),
            DiscreteGaussian::Ziggurat(sampler) => sampler.sample(rng),
        }
    }
}
