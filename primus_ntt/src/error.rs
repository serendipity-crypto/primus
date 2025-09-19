use std::fmt::Debug;

use thiserror::Error;

/// Errors that may occur.
#[derive(Error, Debug)]
pub enum NttError<T> {
    /// Error that occurs when the given modulus has no primitive root with the given degree.
    #[error("There is no primitive root with the degree {degree:?} and the modulus {modulus:?}!")]
    NoPrimitiveRoot {
        /// The degree for the primitive root
        degree: T,
        /// The modulus.
        modulus: T,
    },
}
