mod error;

mod common;

mod binary;
mod ternary;

mod discrete_gaussian;

pub use error::DistrErr;

pub use common::{sample_binary_values, sample_ternary_values};

pub use binary::BinaryDistr;
pub use ternary::TernaryDistr;

#[cfg(target_os = "linux")]
pub use discrete_gaussian::UnixCDTSampler;
pub use discrete_gaussian::{CDTSampler, DiscreteGaussian, DiscreteZiggurat};
