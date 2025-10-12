use primus_integer::{UnsignedInteger, size::Size};
use primus_lattice::Rlwe;
use primus_poly::{ArrayBase, Polynomial};
use primus_reduce::RingContext;
use rand::distr::Distribution;

use crate::{
    LweCiphertext, LweParameters, LweSecretKeyType, MultiMsgLweCiphertext, decode, encode,
};

/// Represents a secret key for the Learning with Errors (LWE) cryptographic scheme.
#[derive(Clone)]
pub struct LweSecretKey<T: UnsignedInteger> {
    pub(crate) key: Vec<T>,
    pub(crate) distr: LweSecretKeyType,
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

        let distr = modulus.uniform_distribution();

        let mut data: Rlwe<Vec<T>> = Rlwe::zero(dimension * 2);

        let (a, b) = data.a_b_mut_slices();

        for (i, o) in a.iter_mut().zip(distr.sample_iter(&mut *csrng)) {
            *i = o;
        }

        Polynomial(ArrayBase(a)).naive_mul_inplace(
            &Polynomial(ArrayBase(self.key.as_ref())),
            &mut Polynomial(ArrayBase(&mut *b)),
            modulus,
        );

        for (&message, bi) in messages.iter().zip(b.iter_mut()) {
            modulus.reduce_add_assign(
                bi,
                encode(message, params.plain_modulus_value, modulus.value()),
            );
        }

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *csrng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        data.extract_first_few_lwe_locally(messages.len(), modulus)
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

        let distr = modulus.uniform_distribution();

        let mut data: Rlwe<Vec<T>> = Rlwe::zero(dimension * 2);

        let (a, b) = data.a_b_mut_slices();

        a.iter_mut()
            .zip(distr.sample_iter(&mut *csrng))
            .for_each(|(i, o): (&mut T, T)| {
                *i = o;
            });

        Polynomial(ArrayBase(a)).naive_mul_inplace(
            &Polynomial(ArrayBase(self.key.as_ref())),
            &mut Polynomial(ArrayBase(&mut *b)),
            modulus,
        );

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *csrng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        data.extract_first_few_lwe_locally(zero_count, modulus)
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
