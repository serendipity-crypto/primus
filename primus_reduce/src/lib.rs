mod common;
mod error;

pub mod lazy_ops;
pub mod ops;

pub use common::{FieldContext, RingContext};
pub use error::ReduceError;

use lazy_ops::*;
use ops::*;

/// Trait for types that represent a modulus.
pub trait Modulus: Copy {
    type ValueT;

    /// Returns the modulus value.
    fn value(self) -> Option<Self::ValueT>;

    /// Returns the modulus value without checking.
    fn value_unchecked(self) -> Self::ValueT;

    /// Returns the value of the modulus minus one.
    fn minus_one(self) -> Self::ValueT;
}
