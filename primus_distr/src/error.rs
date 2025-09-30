use std::fmt::Debug;

use primus_integer::UnsignedInteger;
use thiserror::Error;

/// Errors that may occur.
#[derive(Error, Debug)]
pub enum DistrErr<T: UnsignedInteger> {
    /// Error that occurs when fails to generate the distribution.
    #[error(
        "Fail to generate the desired discrete gaussian distribution.\nmean:{mean}\nstandard deviation:{std_dev}\nmodulus_minus_one:{modulus_minus_one}"
    )]
    DiscreteGaussianErr {
        mean: f64,
        std_dev: f64,
        modulus_minus_one: T,
    },
}
