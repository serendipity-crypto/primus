use primus_integer::{UnsignedInteger, size::Size};
use primus_reduce::RingContext;
use rand::distr::Distribution;
use zeroize::{Zeroize, ZeroizeOnDrop};

use crate::{
    LweCiphertext, LweParameters, LweSecretKeyType, MultiMsgLweCiphertext, decode, encode,
};

/// Represents a secret key for the Learning with Errors (LWE) cryptographic scheme.
#[derive(Clone)]
pub struct LweSecretKey<T: UnsignedInteger> {
    data: Vec<T>,
    distr: LweSecretKeyType,
}

impl<T: UnsignedInteger> Zeroize for LweSecretKey<T> {
    #[inline]
    fn zeroize(&mut self) {
        self.data.zeroize();
    }
}

impl<T: UnsignedInteger> ZeroizeOnDrop for LweSecretKey<T> {}

impl<T: UnsignedInteger> AsRef<[T]> for LweSecretKey<T> {
    #[inline]
    fn as_ref(&self) -> &[T] {
        &self.data
    }
}

impl<T: UnsignedInteger> Size for LweSecretKey<T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.data.byte_count()
    }
}

impl<T: UnsignedInteger> LweSecretKey<T> {
    /// Creates a new [`LweSecretKey<T>`].
    #[inline]
    pub fn new(key: Vec<T>, distr: LweSecretKeyType) -> Self {
        Self { data: key, distr }
    }

    /// Returns the dimension of this [`LweSecretKey<T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.data.len()
    }

    /// Returns the distribution of this [`LweSecretKey<T>`].
    #[inline]
    pub fn distr(&self) -> LweSecretKeyType {
        self.distr
    }

    /// Generates a new [`LweSecretKey<T>`] with random values.
    #[inline]
    pub fn generate<R, M>(params: &LweParameters<T, M>, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: RingContext<T>,
    {
        let distr = params.secret_key_type();
        let key = match distr {
            LweSecretKeyType::Binary => primus_distr::sample_binary_values(params.dimension(), rng),
            LweSecretKeyType::Ternary => primus_distr::sample_ternary_values(
                params.cipher_modulus_minus_one(),
                params.dimension(),
                rng,
            ),
        };
        Self { data: key, distr }
    }

    /// Encrypts message into [`LweCiphertext<T>`].
    #[inline]
    pub fn encrypt<R, M, Msg>(
        &self,
        message: Msg,
        params: &LweParameters<T, M>,
        rng: &mut R,
    ) -> LweCiphertext<T>
    where
        Msg: TryInto<T>,
        R: rand::Rng + rand::CryptoRng,
        M: RingContext<T>,
    {
        let gaussian = params.noise_distribution();
        let modulus = params.cipher_modulus();
        let uniform = params.cipher_modulus_uniform_distr();

        let mut ciphertext = LweCiphertext::generate_random_zero_sample(
            self.as_ref(),
            modulus,
            uniform,
            gaussian,
            rng,
        );
        modulus.reduce_add_assign(
            ciphertext.b_mut(),
            encode(message, params.plain_modulus_value(), modulus.value()),
        );

        ciphertext
    }

    /// Encrypts multiple messages using the secret key.
    #[inline]
    pub fn encrypt_multi_messages<R, M, Msg>(
        &self,
        messages: &[Msg],
        params: &LweParameters<T, M>,
        rng: &mut R,
    ) -> MultiMsgLweCiphertext<T>
    where
        Msg: Copy + TryInto<T>,
        R: rand::Rng + rand::CryptoRng,
        M: RingContext<T>,
    {
        let dimension = params.dimension();
        let gaussian = params.noise_distribution();
        let uniform = params.cipher_modulus_uniform_distr();
        let modulus = params.cipher_modulus();

        let mut data: Vec<T> = vec![T::ZERO; dimension + messages.len()];
        let (a, b) = data.split_at_mut(dimension);

        for (i, o) in a.iter_mut().zip(uniform.sample_iter(&mut *rng)) {
            *i = o;
        }

        b.iter_mut().enumerate().for_each(|(i, bi)| {
            if i != 0 {
                *bi = modulus.reduce_dot_product_iter(
                    a.iter()
                        .skip(dimension - i)
                        .map(|&ai| modulus.reduce_neg(ai))
                        .chain(a.iter().take(dimension - i).copied()),
                    self.data.iter().copied(),
                )
            } else {
                *bi = modulus.reduce_dot_product(&a, self);
            }
        });

        for (bi, &message) in b.iter_mut().zip(messages) {
            modulus.reduce_add_assign(
                bi,
                encode(message, params.plain_modulus_value(), modulus.value()),
            );
        }

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *rng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        MultiMsgLweCiphertext::new(data)
    }

    /// Encrypts multiple zeros using the secret key.
    #[inline]
    pub fn encrypt_multi_zeros<R, Modulus>(
        &self,
        zero_count: usize,
        params: &LweParameters<T, Modulus>,
        rng: &mut R,
    ) -> MultiMsgLweCiphertext<T>
    where
        R: rand::Rng + rand::CryptoRng,
        Modulus: RingContext<T>,
    {
        let dimension = params.dimension();
        let gaussian = params.noise_distribution();
        let uniform = params.cipher_modulus_uniform_distr();
        let modulus = params.cipher_modulus();

        let mut data: Vec<T> = vec![T::ZERO; dimension + zero_count];
        let (a, b) = data.split_at_mut(dimension);

        a.iter_mut()
            .zip(uniform.sample_iter(&mut *rng))
            .for_each(|(i, o): (&mut T, T)| {
                *i = o;
            });

        b.iter_mut().enumerate().for_each(|(i, bi)| {
            if i != 0 {
                *bi = modulus.reduce_dot_product_iter(
                    a.iter()
                        .skip(dimension - i)
                        .map(|&ai| modulus.reduce_neg(ai))
                        .chain(a.iter().take(dimension - i).copied()),
                    self.data.iter().copied(),
                )
            } else {
                *bi = modulus.reduce_dot_product(&a, self);
            }
        });

        for (bi, ei) in b.iter_mut().zip(gaussian.sample_iter(&mut *rng)) {
            modulus.reduce_add_assign(bi, ei);
        }

        MultiMsgLweCiphertext::new(data)
    }

    /// Decrypts the [`LweCiphertext<T>`] back to message.
    #[inline]
    pub fn decrypt<M, Msg>(
        &self,
        cipher_text: &LweCiphertext<T>,
        params: &LweParameters<T, M>,
    ) -> Msg
    where
        Msg: TryFrom<T>,
        M: RingContext<T>,
    {
        let modulus = params.cipher_modulus();

        let (a, b) = cipher_text.a_b();

        let a_mul_s = modulus.reduce_dot_product(a, self);
        let plaintext = modulus.reduce_sub(b, a_mul_s);

        decode(plaintext, params.plain_modulus_value(), modulus.value())
    }

    /// Decrypts the [`LweCiphertext<T>`] back to message.
    #[inline]
    pub fn decrypt_with_noise<M, Msg>(
        &self,
        cipher_text: &LweCiphertext<T>,
        params: &LweParameters<T, M>,
    ) -> (Msg, T)
    where
        Msg: Copy + TryFrom<T> + TryInto<T>,
        M: RingContext<T>,
    {
        let modulus = params.cipher_modulus();

        let (a, b) = cipher_text.a_b();

        let a_mul_s = modulus.reduce_dot_product(a, self);
        let plaintext = modulus.reduce_sub(b, a_mul_s);

        let t = params.plain_modulus_value();
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

    /// Decrypts the [`MultiMsgLweCiphertext<T>`] back to message.
    #[inline]
    pub fn decrypt_multi_messages<M, Msg>(
        &self,
        cipher_text: &MultiMsgLweCiphertext<T>,
        params: &LweParameters<T, M>,
    ) -> Vec<Msg>
    where
        Msg: TryFrom<T>,
        M: RingContext<T>,
    {
        let modulus = params.cipher_modulus();
        let dimension = params.dimension();
        let (a, b) = cipher_text.a_b(dimension);

        b.iter()
            .enumerate()
            .map(|(i, &b)| {
                let a_mul_s = modulus.reduce_dot_product_iter(
                    a.iter()
                        .skip(dimension - i)
                        .map(|&v| modulus.reduce_neg(v))
                        .chain(a.iter().take(dimension - i).copied()),
                    self.data.iter().copied(),
                );
                let plaintext = modulus.reduce_sub(b, a_mul_s);

                decode(plaintext, params.plain_modulus_value(), modulus.value())
            })
            .collect()
    }
}
