use std::ops::Deref;

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
pub struct RlweSecretKey<T>
where
    T: UnsignedInteger,
{
    key: PolynomialOwned<T>,
    distr: RingSecretKeyType,
}

impl<T> Deref for RlweSecretKey<T>
where
    T: UnsignedInteger,
{
    type Target = PolynomialOwned<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T> RlweSecretKey<T>
where
    T: UnsignedInteger,
{
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
    pub fn generate<R, M>(params: &RlweParameters<T, M>, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let distr = params.secret_key_type();
        let poly_length = params.poly_length();
        let modulus_minus_one = params.cipher_modulus_minus_one();

        let key = match distr {
            RingSecretKeyType::Binary => Polynomial::random_binary(poly_length, rng),
            RingSecretKeyType::Ternary => {
                Polynomial::random_ternary(modulus_minus_one, poly_length, rng)
            }
            RingSecretKeyType::Gaussian => {
                Polynomial::random_gaussian(poly_length, params.noise_distribution(), rng)
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
pub struct NttRlweSecretKey<T>
where
    T: UnsignedInteger,
{
    key: NttPolynomialOwned<T>,
    distr: RingSecretKeyType,
}

impl<T> Deref for NttRlweSecretKey<T>
where
    T: UnsignedInteger,
{
    type Target = NttPolynomialOwned<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T> NttRlweSecretKey<T>
where
    T: UnsignedInteger,
{
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
        let key = secret_key.key.clone();
        let key = ntt_table.transform_inplace(key);
        Self {
            key,
            distr: secret_key.distr,
        }
    }

    /// Performs `b-as`.
    pub fn phase_inplace<Table, M, A>(
        &self,
        cipher: &NttRlweCiphertext<A>,
        result: &mut PolynomialOwned<T>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        A: RawData<Elem = T> + Data,
    {
        let (a, b) = cipher.a_b_slices();

        NttPolynomial(ArrayBase(a)).mul_inplace(
            &self.key,
            &mut NttPolynomial(ArrayBase(result.as_mut())),
            modulus,
        );
        NttPolynomial(ArrayBase(b))
            .sub_to_right(&mut NttPolynomial(ArrayBase(result.as_mut())), modulus);

        ntt_table.inverse_transform_slice(result.as_mut())
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
        TruncatedRlwe::generate_random_zero_sample(
            zero_count,
            &self.key,
            params.cipher_modulus_uniform_distr(),
            params.noise_distribution(),
            ntt_table,
            params.cipher_modulus(),
            rng,
        )
    }

    /// Decrypts the [`TruncatedRlwe`] back to message.
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
        let modulus_value = modulus.value();

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
                    params.plain_modulus_value(),
                    modulus_value,
                )
            })
            .collect()
    }

    /// Decrypts the [`TruncatedRlwe`] back to message.
    #[inline]
    pub fn phase_multi_messages<M, Table>(
        &self,
        cipher: &TruncatedRlwe<Vec<T>>,
        modulus: M,
        ntt_table: &Table,
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
