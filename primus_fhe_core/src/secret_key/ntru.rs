use std::ops::Deref;

use primus_integer::{Data, RawData, UnsignedInteger};
use primus_ntt::NttTable;
use primus_poly::{NttPolynomial, NttPolynomialOwned, PolynomialOwned};
use primus_reduce::FieldContext;
use zeroize::{Zeroize, ZeroizeOnDrop};

use crate::{NttNtruCiphertext, RingSecretKeyType};

/// Represents a secret key for the NTRU cryptographic scheme.
#[derive(Clone)]
pub struct NtruSecretKey<T: UnsignedInteger> {
    key: PolynomialOwned<T>,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Zeroize for NtruSecretKey<T> {
    #[inline]
    fn zeroize(&mut self) {
        self.key.0.zeroize();
    }
}

impl<T: UnsignedInteger> ZeroizeOnDrop for NtruSecretKey<T> {}

impl<T: UnsignedInteger> Deref for NtruSecretKey<T> {
    type Target = PolynomialOwned<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T: UnsignedInteger> NtruSecretKey<T> {
    /// Creates a new [`NtruSecretKey<T>`].
    pub fn new(key: PolynomialOwned<T>, distr: RingSecretKeyType) -> Self {
        Self { key, distr }
    }

    /// Returns the distribution of this [`NtruSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }
}

/// Represents a secret key for the (NTT) NTRU cryptographic scheme.
#[derive(Clone)]
pub struct NttNtruSecretKey<T: UnsignedInteger> {
    key: NttPolynomialOwned<T>,
    inv_key: NttPolynomialOwned<T>,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Zeroize for NttNtruSecretKey<T> {
    #[inline]
    fn zeroize(&mut self) {
        self.key.0.zeroize();
        self.inv_key.0.zeroize();
    }
}

impl<T: UnsignedInteger> ZeroizeOnDrop for NttNtruSecretKey<T> {}

impl<T: UnsignedInteger> NttNtruSecretKey<T> {
    /// Creates a new [`NttNtruSecretKey<T>`].
    pub fn new(
        key: NttPolynomialOwned<T>,
        inv_key: NttPolynomialOwned<T>,
        distr: RingSecretKeyType,
    ) -> Self {
        Self {
            key,
            inv_key,
            distr,
        }
    }

    /// Returns the distribution of this [`NttNtruSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    /// Performs `h*f`.
    pub fn phase_inplace<Table, M, A>(
        &self,
        cipher: &NttNtruCiphertext<A>,
        result: &mut PolynomialOwned<T>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T>,
        A: RawData<Elem = T> + Data,
    {
        let h = cipher.as_ref();
        let mut temp = NttPolynomial(result.as_mut());
        NttPolynomial(h).mul_inplace(&self.key, &mut temp, modulus);
        ntt_table.inverse_transform_slice(result.as_mut())
    }
}
