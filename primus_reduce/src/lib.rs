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
//!
//! # Implementing [`RingContext`] / [`FieldContext`]
//!
//! Both are *marker* traits with blanket impls: implement every listed
//! `Reduce*` (and `LazyReduce*`) trait for your modulus type and the
//! corresponding context trait is granted automatically.

use core::fmt::Debug;

mod common;
mod error;

pub mod lazy_ops;
pub mod ops;

pub use common::{FieldContext, RingContext};
pub use error::ReduceError;
pub use lazy_ops::*;
pub use ops::*;

use num_traits::ConstZero;
use primus_integer::Integer;
use rand::distr::Uniform;

/// Trait for types that represent a modulus.
pub trait Modulus: Copy + Debug + Send + Sync {
    type ValueT: Integer;

    /// Returns the modulus value, or `None` when the modulus is implicit
    /// (e.g. a native power-of-two modulus where the value is `2^BITS` and
    /// cannot be represented in `ValueT`).
    #[must_use]
    fn value(self) -> Option<Self::ValueT>;

    /// Returns the modulus value without checking that it fits in `ValueT`.
    ///
    /// # Returns
    ///
    /// For an explicit modulus, returns the same value as [`value`](Self::value)
    /// would unwrap. For an implicit modulus (e.g. native `2^BITS`), the
    /// result is implementation-defined (typically `0` from wrapping); use
    /// [`value`](Self::value) when correctness depends on the distinction.
    #[must_use]
    fn value_unchecked(self) -> Self::ValueT;

    /// Returns the value of the modulus minus one.
    ///
    /// Well-defined for both explicit and implicit moduli: for the implicit
    /// native power-of-two case this is `T::MAX`.
    #[must_use]
    fn minus_one(self) -> Self::ValueT;

    /// Returns a [Uniform] distribution over the values of [Modulus].
    ///
    /// # Panics
    ///
    /// Never panics for unsigned [`Modulus`] impls: the constructed range
    /// `[0, minus_one()]` is always non-empty.
    #[must_use]
    #[inline]
    fn uniform_distribution(self) -> Uniform<Self::ValueT> {
        Uniform::new_inclusive(<Self::ValueT as ConstZero>::ZERO, self.minus_one())
            .expect("uniform_distribution: invalid modulus range")
    }
}
