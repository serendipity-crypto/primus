use core::fmt::Debug;

use thiserror::Error;

/// Errors that may occur.
#[derive(Error, Debug)]
pub enum ReduceError<T> {
    /// Error that occurs when the given value has no inverse element with the given modulus.
    #[error("Value {value:?} has no inverse element with the modulus {modulus:?}!")]
    NoInverse {
        /// The value being inverted.
        value: T,
        /// The modulus.
        modulus: T,
    },
}
