mod error;

mod common;

mod binary;
mod ternary;

mod discrete_gaussian;
mod signed_discrete_gaussian;

pub use error::DistrErr;

pub use common::*;

pub use binary::BinaryDistr;
pub use ternary::TernaryDistr;

#[cfg(all(target_os = "linux", feature = "high_precision"))]
pub use discrete_gaussian::UnixCDTSampler;
pub use discrete_gaussian::{CDTSampler, DiscreteGaussian, DiscreteZiggurat};
pub use signed_discrete_gaussian::SignedDiscreteGaussian;
