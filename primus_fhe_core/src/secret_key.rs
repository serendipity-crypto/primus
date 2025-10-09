use std::ops::Deref;

use primus_distr::DiscreteGaussian;
use primus_integer::{UnsignedInteger, size::Size};
use primus_lattice::Rlwe;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{NttPolynomial, Polynomial};
use primus_reduce::{FieldContext, RingContext};
use rand::distr::{Distribution, Uniform};

use crate::{
    LweCiphertext, LweParameters, MultiMsgLweCiphertext, NttRlweCiphertext, decode, encode,
};

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
#[derive(Clone)]
pub struct LweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    distr: LweSecretKeyType,
}

impl<T: UnsignedInteger> AsRef<[T]> for LweSecretKey<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        &self.key
    }
}

impl<T: UnsignedInteger> Size for LweSecretKey<T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.key.byte_count()
    }
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

    /// Generates a new `LweSecretKey` with random values.
    #[inline]
    pub fn generate<R, M>(params: &LweParameters<T, M>, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: RingContext<T>,
    {
        let distr = params.secret_key_type;
        let key = match distr {
            LweSecretKeyType::Binary => primus_distr::sample_binary_values(params.dimension, rng),
            LweSecretKeyType::Ternary => primus_distr::sample_ternary_values(
                params.cipher_modulus_minus_one,
                params.dimension,
                rng,
            ),
        };
        Self { key, distr }
    }

    /// Encrypts message into [`LweCiphertext<C>`].
    #[inline]
    pub fn encrypt<Msg, R, Modulus>(
        &self,
        message: Msg,
        params: &LweParameters<T, Modulus>,
        rng: &mut R,
    ) -> LweCiphertext<T>
    where
        Msg: TryInto<T>,
        R: rand::Rng + rand::CryptoRng,
        Modulus: RingContext<T>,
    {
        let gaussian = params.noise_distribution();
        let modulus = params.cipher_modulus;

        let mut ciphertext =
            LweCiphertext::generate_random_zero_sample(self.as_ref(), modulus, &gaussian, rng);
        modulus.reduce_add_assign(
            ciphertext.b_mut(),
            encode(message, params.plain_modulus_value, modulus.value()),
        );

        ciphertext
    }

    /// Encrypts multiple messages using the secret key.
    ///
    /// # Arguments
    ///
    /// * `messages` - A slice of messages to be encrypted.
    /// * `params` - The parameters for the LWE scheme.
    /// * `csrng` - A mutable reference to a random number generator.
    ///
    /// # Returns
    ///
    /// A `CmLweCiphertext` containing the encrypted messages.
    #[inline]
    pub fn encrypt_multi_messages<Msg, R, Modulus>(
        &self,
        messages: &[Msg],
        params: &LweParameters<T, Modulus>,
        csrng: &mut R,
    ) -> MultiMsgLweCiphertext<T>
    where
        Msg: Copy + TryInto<T>,
        R: rand::Rng + rand::CryptoRng,
        Modulus: RingContext<T>,
    {
        let dimension = params.dimension;
        let gaussian = params.noise_distribution();
        let modulus = params.cipher_modulus;

        let distr =
            <Uniform<T>>::new_inclusive(T::ZERO, params.cipher_modulus_minus_one()).unwrap();

        let mut a = Polynomial::zero(dimension);
        let mut b = Polynomial::zero(dimension);

        for (i, o) in a.iter_mut().zip(distr.sample_iter(&mut *csrng)) {
            *i = o;
        }

        a.naive_mul_inplace(&self.key, &mut b, modulus);

        for (&message, bi) in messages.iter().zip(b.iter_mut()) {
            modulus.reduce_add_assign(
                bi,
                encode(message, params.plain_modulus_value, modulus.value()),
            );
        }

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *csrng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        Rlwe::new(a, b).extract_first_few_lwe_locally(messages.len(), modulus)
    }

    /// Encrypts multiple zeros using the secret key.
    ///
    /// # Arguments
    ///
    /// * `zero_count` - The count of zeros to be encrypted.
    /// * `params` - The parameters for the LWE scheme.
    /// * `csrng` - A mutable reference to a random number generator.
    ///
    /// # Returns
    ///
    /// A `MultiMsgLweCiphertext` containing the encrypted messages.
    #[inline]
    pub fn encrypt_multi_zeros<R, Modulus>(
        &self,
        zero_count: usize,
        params: &LweParameters<T, Modulus>,
        csrng: &mut R,
    ) -> MultiMsgLweCiphertext<T>
    where
        R: rand::Rng + rand::CryptoRng,
        Modulus: RingContext<T>,
    {
        let dimension = params.dimension;
        let gaussian = params.noise_distribution();
        let modulus = params.cipher_modulus;

        let distr =
            <Uniform<T>>::new_inclusive(T::ZERO, params.cipher_modulus_minus_one()).unwrap();

        let mut a = Polynomial::zero(dimension);
        let mut b = Polynomial::zero(dimension);

        a.iter_mut()
            .zip(distr.sample_iter(&mut *csrng))
            .for_each(|(i, o): (&mut T, T)| {
                *i = o;
            });

        a.naive_mul_inplace(&self.key, &mut b, modulus);

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *csrng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        Rlwe::new(a, b).extract_first_few_lwe_locally(zero_count, modulus)
    }

    /// Decrypts the [`LweCiphertext`] back to message.
    #[inline]
    pub fn decrypt<Msg, Modulus>(
        &self,
        cipher_text: &LweCiphertext<T>,
        params: &LweParameters<T, Modulus>,
    ) -> Msg
    where
        Msg: TryFrom<T>,
        Modulus: RingContext<T>,
    {
        let modulus = params.cipher_modulus;

        let a_mul_s = modulus.reduce_dot_product(cipher_text.a(), self);
        let plaintext = modulus.reduce_sub(cipher_text.b(), a_mul_s);

        decode(plaintext, params.plain_modulus_value, modulus.value())
    }

    /// Decrypts the [`LweCiphertext`] back to message.
    #[inline]
    pub fn decrypt_with_noise<Msg, Modulus>(
        &self,
        cipher_text: &LweCiphertext<T>,
        params: &LweParameters<T, Modulus>,
    ) -> (Msg, T)
    where
        Msg: Copy + TryFrom<T> + TryInto<T>,
        Modulus: RingContext<T>,
    {
        let modulus = params.cipher_modulus;
        let a_mul_s = modulus.reduce_dot_product(cipher_text.a(), self);
        let plaintext = modulus.reduce_sub(cipher_text.b(), a_mul_s);

        let t = params.plain_modulus_value;
        let q = modulus.value();
        let message = decode(plaintext, t, q);
        let fresh = encode(message, t, q);

        (
            message,
            modulus
                .reduce_sub(plaintext, fresh)
                .min(modulus.reduce_sub(fresh, plaintext)),
        )
    }

    /// Decrypts the [`LweCiphertext`] back to message.
    #[inline]
    pub fn decrypt_multi_messages<Msg, Modulus>(
        &self,
        cipher_text: &MultiMsgLweCiphertext<T>,
        params: &LweParameters<T, Modulus>,
    ) -> Vec<Msg>
    where
        Msg: TryFrom<T>,
        Modulus: RingContext<T>,
    {
        let modulus = params.cipher_modulus;
        let dimension = cipher_text.a().len();

        cipher_text
            .b()
            .iter()
            .enumerate()
            .map(|(i, &b)| {
                let a_mul_s = modulus.reduce_dot_product_iter(
                    cipher_text.a()[dimension - i..]
                        .iter()
                        .map(|&v| modulus.reduce_neg(v))
                        .chain(cipher_text.a()[..dimension - i].iter().copied()),
                    self.key.iter().copied(),
                );
                let plaintext = modulus.reduce_sub(b, a_mul_s);

                decode(plaintext, params.plain_modulus_value, modulus.value())
            })
            .collect()
    }
}

/// Represents a secret key for the Ring Learning with Errors (RLWE) cryptographic scheme.
#[derive(Clone)]
pub struct RlweSecretKey<T: UnsignedInteger> {
    key: Polynomial<T>,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Deref for RlweSecretKey<T> {
    type Target = Polynomial<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T: UnsignedInteger> RlweSecretKey<T> {
    /// Creates a new [`RlweSecretKey<T>`].
    #[inline]
    pub fn new(key: Polynomial<T>, distr: RingSecretKeyType) -> Self {
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
        let distr = match lwe_secret_key.distr {
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
    key: NttPolynomial<T>,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> Deref for NttRlweSecretKey<T> {
    type Target = NttPolynomial<T>;

    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.key
    }
}

impl<T: UnsignedInteger> NttRlweSecretKey<T> {
    /// Creates a new [`NttRlweSecretKey<T>`].
    #[inline]
    pub fn new(key: NttPolynomial<T>, distr: RingSecretKeyType) -> Self {
        Self { key, distr }
    }

    /// Returns the distribution type of the secret key.
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    /// Creates a new `NttRlweSecretKey` from a coefficient secret key.
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &RlweSecretKey<T>, ntt_table: &Table) -> Self
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        Self {
            key: ntt_table.transform(&secret_key.key),
            distr: secret_key.distr,
        }
    }

    /// Performs `b-as`.
    pub fn phase<Table, M>(
        &self,
        cipher: NttRlweCiphertext<T>,
        ntt_table: &Table,
        modulus: M,
    ) -> Polynomial<T>
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a_ntt, b_ntt) = cipher.into_inner();
        let a_mul_s_ntt = a_ntt.mul(&self.key, modulus);
        let dec_poly_ntt = b_ntt.sub(&a_mul_s_ntt, modulus);
        ntt_table.inverse_transform_inplace(dec_poly_ntt)
    }

    // /// Encrypts multiple zeros using the secret key.
    // pub fn encrypt_multi_zeros<R>(
    //     &self,
    //     zero_count: usize,
    //     params: &RlweParameters<T>,
    //     ntt_table: &<F as NttField>::Table,
    //     rng: &mut R,
    // ) -> SparseRlwe<F>
    // where
    //     R: Rng + CryptoRng,
    // {
    //     SparseRlwe::generate_random_zero_sample(
    //         zero_count,
    //         &self.key,
    //         params.noise_distribution(),
    //         ntt_table,
    //         rng,
    //     )
    // }
}

/// Represents a secret key for the Ring Learning with Errors (RLWE) cryptographic scheme.
#[derive(Clone)]
pub struct GlweSecretKey<T: UnsignedInteger> {
    key: Vec<Polynomial<T>>,
    distr: RingSecretKeyType,
}
