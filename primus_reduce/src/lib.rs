//! Traits for modular reduction operations.
//!
//! This crate defines the algebraic interface between modulus types
//! (implemented in `primus_modulus`) and values.  Each arithmetic
//! operation — reduce, add, sub, mul, square, neg, exp, dot-product,
//! inverse, division — has its own fine-grained trait so that modulus
//! types only need to implement the subset they actually support.
//!
//! The two marker supertraits [`RingContext`] and [`FieldContext`]
//! aggregate the full ring / field operation sets respectively.

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

    /// Returns the modulus value, or `None` when the modulus is implicit
    /// (e.g. a native power-of-two modulus where the value is `2^BITS` and
    /// cannot be represented in `ValueT`).
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
