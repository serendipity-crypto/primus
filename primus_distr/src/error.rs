use std::fmt::Debug;

use primus_integer::Integer;
use thiserror::Error;

/// Errors that may occur.
#[derive(Error, Debug)]
pub enum DistrErr<T: Integer> {
    /// Error that occurs when fails to generate the distribution.
    #[error(
        "Fail to generate the desired discrete gaussian distribution.\nstandard deviation:{std_dev}\nmodulus_minus_one:{modulus_minus_one}"
    )]
    DiscreteGaussianErr { std_dev: f64, modulus_minus_one: T },
}
