mod common;
mod error;

pub mod lazy_ops;
pub mod ops;

pub use common::{FieldContext, RingContext};
pub use error::ReduceError;

use lazy_ops::*;
use num_traits::ConstZero;
use ops::*;
use primus_integer::Integer;
use rand::distr::Uniform;

/// Trait for types that represent a modulus.
pub trait Modulus: Copy {
    type ValueT: Integer;

    /// Returns the modulus value.
    fn value(self) -> Option<Self::ValueT>;

    /// Returns the modulus value without checking.
    fn value_unchecked(self) -> Self::ValueT;

    /// Returns the value of the modulus minus one.
    fn minus_one(self) -> Self::ValueT;

    /// Returns a [Uniform] distribution over the values of [Modulus].
    #[must_use]
    #[inline]
    fn uniform_distribution(self) -> Uniform<Self::ValueT> {
        Uniform::new_inclusive(<Self::ValueT as ConstZero>::ZERO, self.minus_one()).unwrap()
    }
}
