//! Value-side modular arithmetic traits.
//!
//! This crate is the mirror of [`primus_reduce`]: where `primus_reduce`
//! attaches operations to the *modulus* (`modulus.reduce_add(a, b)`),
//! this crate attaches them to the *value* (`a.add_modulo(b, modulus)`).
//! Each trait is implemented via a blanket impl that delegates to the
//! corresponding [`primus_reduce`] trait, simply reversing the call order.
//!
//! The naming convention is `XxxModulo` (value-side) ↔ `ReduceXxx` (modulus-side).

mod error;

pub mod lazy_ops;
pub mod ops;

pub use error::ModuloError;

pub use lazy_ops::*;
pub use ops::*;
