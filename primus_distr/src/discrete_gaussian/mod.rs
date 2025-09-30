use primus_integer::UnsignedInteger;
use rand::{Rng, distr::Distribution};

mod cdt;
#[cfg(target_os = "linux")]
mod unix_cdt;
mod ziggurat;

pub use cdt::CDTSampler;
#[cfg(target_os = "linux")]
pub use unix_cdt::UnixCDTSampler;
pub use ziggurat::DiscreteZiggurat;

use crate::DistrErr;

/// The gaussian distribution `N(mean, std_dev**2)`.
#[derive(Clone)]
pub enum DiscreteGaussian<T: UnsignedInteger> {
    /// CDTSampler
    Cdt(super::CDTSampler<T>),
    /// UnixCDTSampler
    #[cfg(target_os = "linux")]
    Unix(super::UnixCDTSampler<T>),
    /// DiscreteZiggurat
    Ziggurat(super::DiscreteZiggurat<T>),
}

impl<T: UnsignedInteger> DiscreteGaussian<T> {
    /// Construct, from mean and standard deviation
    ///
    /// Parameters:
    ///
    /// -   mean (`μ`, unrestricted)
    /// -   standard deviation (`σ`, must be finite)
    #[inline]
    pub fn new(
        _mean: f64,
        std_dev: f64,
        modulus_minus_one: T,
    ) -> Result<DiscreteGaussian<T>, DistrErr<T>> {
        if std_dev < 3.0 {
            #[cfg(target_os = "linux")]
            {
                Ok(DiscreteGaussian::Unix(super::UnixCDTSampler::new(
                    std_dev,
                    6.0,
                    modulus_minus_one,
                )))
            }

            #[cfg(not(target_os = "linux"))]
            Ok(DiscreteGaussian::Cdt(super::CDTSampler::new(
                std_dev,
                6.0,
                modulus_minus_one,
            )))
        } else {
            Ok(DiscreteGaussian::Ziggurat(super::DiscreteZiggurat::new(
                std_dev,
                6.0,
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
        mean: f64,
        std_dev: f64,
        max_std_dev: f64,
        modulus_minus_one: T,
    ) -> Result<DiscreteGaussian<T>, DistrErr<T>> {
        if max_std_dev <= std_dev || std_dev < 0.7 {
            return Err(DistrErr::DiscreteGaussianErr {
                mean,
                std_dev,
                modulus_minus_one,
            });
        }
        unimplemented!()
    }
}

impl<T: UnsignedInteger> Distribution<T> for DiscreteGaussian<T> {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> T {
        match self {
            DiscreteGaussian::Cdt(cdtsampler) => cdtsampler.sample(rng),
            #[cfg(target_os = "linux")]
            DiscreteGaussian::Unix(sampler) => sampler.sample(rng),
            DiscreteGaussian::Ziggurat(discrete_ziggurat) => discrete_ziggurat.sample(rng),
        }
    }
}
