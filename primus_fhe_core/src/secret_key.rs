// mod glwe;
// mod lwe;
// mod rlwe;

use primus_integer::UnsignedInteger;

/// The distribution type of the LWE Secret Key.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum LweSecretKeyType {
    /// Binary SecretKey Distribution.
    Binary,
    /// Ternary SecretKey Distribution.
    #[default]
    Ternary,
}

/// The distribution type of the Ring Secret Key.
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum RingSecretKeyType {
    /// Binary SecretKey Distribution.
    Binary,
    /// Ternary SecretKey Distribution.
    #[default]
    Ternary,
    /// Gaussian SecretKey Distribution.
    Gaussian,
}

/// Represents a secret key for the Learning with Errors (LWE) cryptographic scheme.
///
/// # Type Parameters
///
/// * `T` - An unsigned integer type that represents the coefficients of the secret key.
#[derive(Clone)]
pub struct LweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    distr: LweSecretKeyType,
}

impl<T: UnsignedInteger> LweSecretKey<T> {
    /// Creates a new [`LweSecretKey<T>`].
    #[inline]
    pub fn new(key: Vec<T>, distr: LweSecretKeyType) -> Self {
        Self { key, distr }
    }

    /// Returns the dimension of the secret key.
    #[inline]
    pub fn dimension(&self) -> usize {
        self.key.len()
    }

    /// Returns the distribution of this [`LweSecretKey<T>`].
    #[inline]
    pub fn distr(&self) -> LweSecretKeyType {
        self.distr
    }
}

impl<T: UnsignedInteger> AsRef<[T]> for LweSecretKey<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        &self.key
    }
}
