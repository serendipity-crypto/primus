use std::fmt::Debug;

use thiserror::Error;

/// Errors that may occur.
#[derive(Error, Debug)]
pub enum RNSError {
    /// Error that occurs when the given moduli are not coprime.
    #[error("Some moduli are not coprime!")]
    CoPrimeError,
}
