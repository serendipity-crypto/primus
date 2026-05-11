use core::fmt::Debug;

use thiserror::Error;

/// Error type for value-side modulo operations that can fail.
///
/// Wraps [`ReduceError`](primus_reduce::ReduceError); currently only covers
/// the case where a multiplicative inverse does not exist.
#[derive(Error, Debug)]
pub enum ModuloError<T> {
    /// Error that occurs when the given value has no inverse element with the given modulus.
    #[error("Value {value:?} has no inverse element with the modulus {modulus:?}!")]
    NoInverse {
        /// The value being inverted.
        value: T,
        /// The modulus.
        modulus: T,
    },
}
