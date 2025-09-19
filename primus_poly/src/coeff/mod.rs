use integer::UnsignedInteger;
use serde::{Deserialize, Serialize};

/// Represents a polynomial where coefficients are elements of a specified numeric `T`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Polynomial<T> {
    poly: Vec<T>,
}

impl<T> Default for Polynomial<T> {
    #[inline]
    fn default() -> Self {
        Self { poly: Vec::new() }
    }
}
