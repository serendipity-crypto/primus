use std::ops::Deref;

use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_lattice::rlwe::TruncatedRlwe;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{
    ArrayBase, Data, NttPolynomial, NttPolynomialOwned, Polynomial, PolynomialOwned, RawData,
};
use primus_reduce::FieldContext;

use crate::{NttRlweCiphertext, RlweParameters, decode};

use super::{LweSecretKey, LweSecretKeyType, RingSecretKeyType};

/// Represents a secret key for the Ring Learning with Errors (RLWE) cryptographic scheme.
#[derive(Clone)]
pub struct RlweSecretKey<T: UnsignedInteger> {
    pub(crate) key: PolynomialOwned<T>,
    pub(crate) distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Deref for RlweSecretKey<T> {
    type Target = PolynomialOwned<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T: UnsignedInteger> RlweSecretKey<T> {
    /// Creates a new [`RlweSecretKey<T>`].
    #[inline]
    pub fn new(key: PolynomialOwned<T>, distr: RingSecretKeyType) -> Self {
        Self { key, distr }
    }

    /// Returns the distribution type of the secret key.
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    #[inline]
    pub fn generate<R>(
        secret_key_type: RingSecretKeyType,
        poly_length: usize,
        modulus_minus_one: T,
        gaussian: Option<DiscreteGaussian<T>>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        let distr = secret_key_type;
        let key = match distr {
            RingSecretKeyType::Binary => Polynomial::random_binary(poly_length, rng),
            RingSecretKeyType::Ternary => {
                Polynomial::random_ternary(modulus_minus_one, poly_length, rng)
            }
            RingSecretKeyType::Gaussian => {
                Polynomial::random_gaussian(poly_length, &gaussian.unwrap(), rng)
            }
        };

        Self { key, distr }
    }

    #[inline]
    pub fn from_lwe_secret_key<C: UnsignedInteger>(
        lwe_secret_key: &LweSecretKey<C>,
        modulus_minus_one: T,
    ) -> Self {
        let convert = |v: &C| {
            if v.is_zero() {
                T::ZERO
            } else if v.is_one() {
                T::ONE
            } else {
                modulus_minus_one
            }
        };
        let distr = match lwe_secret_key.distr() {
            LweSecretKeyType::Binary => RingSecretKeyType::Binary,
            LweSecretKeyType::Ternary => RingSecretKeyType::Ternary,
        };

        RlweSecretKey {
            key: Polynomial::new(lwe_secret_key.as_ref().iter().map(convert).collect()),
            distr,
        }
    }
}

/// Represents a secret key for the (NTT) Ring Learning with Errors (RLWE) cryptographic scheme.
#[derive(Clone)]
pub struct NttRlweSecretKey<T: UnsignedInteger> {
    key: NttPolynomialOwned<T>,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Deref for NttRlweSecretKey<T> {
    type Target = NttPolynomialOwned<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T: UnsignedInteger> NttRlweSecretKey<T> {
    /// Creates a new [`NttRlweSecretKey<T>`].
    #[inline]
    pub fn new(key: NttPolynomialOwned<T>, distr: RingSecretKeyType) -> Self {
        Self { key, distr }
    }

    /// Returns the distribution type of the secret key.
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    /// Creates a new [`NttRlweSecretKey`] from a coefficient secret key.
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &RlweSecretKey<T>, ntt_table: &Table) -> Self
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let mut data = secret_key.key.0.clone();
        ntt_table.transform_slice(data.as_mut());
        Self {
            key: NttPolynomial(data),
            distr: secret_key.distr,
        }
    }

    /// Performs `b-as`.
    pub fn phase_inplace<Table, M, S>(
        &self,
        cipher: &NttRlweCiphertext<S>,
        result: &mut PolynomialOwned<T>,
        ntt_table: &Table,
        modulus: M,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        S: RawData<Elem = T> + Data,
    {
        let (a, b) = cipher.a_b_slices();

        NttPolynomial(ArrayBase(a)).mul_inplace(
            &self.key,
            &mut NttPolynomial(ArrayBase(result.0.as_mut())),
            modulus,
        );
        NttPolynomial(ArrayBase(b))
            .sub_to_right(&mut NttPolynomial(ArrayBase(result.0.as_mut())), modulus);

        ntt_table.inverse_transform_slice(result.0.as_mut())
    }

    /// Encrypts multiple zeros using the secret key.
    pub fn encrypt_multi_zeros<R, M, Table>(
        &self,
        zero_count: usize,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) -> TruncatedRlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let modulus = params.cipher_modulus();
        TruncatedRlwe::generate_random_zero_sample(
            zero_count,
            &self.key,
            params.noise_distribution(),
            ntt_table,
            modulus,
            rng,
        )
    }

    /// Decrypts the [`LweCiphertext`] back to message.
    #[inline]
    pub fn decrypt_multi_messages<Msg, M, Table>(
        &self,
        cipher: &TruncatedRlwe<Vec<T>>,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
    ) -> Vec<Msg>
    where
        Msg: TryFrom<T>,
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = params.poly_length();
        let modulus = params.cipher_modulus();
        let modulus_value = modulus.value_unchecked();

        let (a, b) = cipher.a_b_slices(poly_length);

        let mut a = a.to_vec();
        ntt_table.transform_slice(&mut a);
        NttPolynomial(ArrayBase(a.as_mut())).mul_assign(&self.key, modulus);
        ntt_table.inverse_transform_slice(&mut a);

        b.iter()
            .zip(a)
            .map(|(x, y)| {
                decode(
                    modulus.reduce_sub(*x, y),
                    {
                        let this = &params;
                        this.plain_modulus_value()
                    },
                    Some(modulus_value),
                )
            })
            .collect()
    }

    /// Decrypts the [`LweCiphertext`] back to message.
    #[inline]
    pub fn phase_multi_messages<M, Table>(
        &self,
        cipher: &TruncatedRlwe<Vec<T>>,
        ntt_table: &Table,
        modulus: M,
    ) -> Vec<T>
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();

        let (a, b) = cipher.a_b_slices(poly_length);

        let mut a = a.to_vec();
        ntt_table.transform_slice(&mut a);
        NttPolynomial(ArrayBase(a.as_mut())).mul_assign(&self.key, modulus);
        ntt_table.inverse_transform_slice(&mut a);

        b.iter()
            .zip(a)
            .map(|(x, y)| modulus.reduce_sub(*x, y))
            .collect()
    }
}
