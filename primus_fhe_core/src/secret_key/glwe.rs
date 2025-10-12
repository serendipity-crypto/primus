use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, NttPolynomial, PolynomialOwned, RawData};
use primus_reduce::FieldContext;

use crate::NttGlweCiphertext;

use super::RingSecretKeyType;

/// Represents a secret key for the Module Learning with Errors (MLWE) cryptographic scheme.
#[derive(Clone)]
pub struct GlweSecretKey<T: UnsignedInteger> {
    pub(crate) key: Vec<T>,
    poly_length: usize,
    dimension: usize,
    pub(crate) distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> GlweSecretKey<T> {
    /// Creates a new [`GlweSecretKey<T>`].
    #[inline]
    pub fn new(
        key: Vec<T>,
        poly_length: usize,
        dimension: usize,
        distr: RingSecretKeyType,
    ) -> Self {
        Self {
            key,
            poly_length,
            dimension,
            distr,
        }
    }

    /// Returns the poly length of this [`GlweSecretKey<T>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the distr of this [`GlweSecretKey<T>`].
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    #[inline]
    pub fn generate<R>(
        secret_key_type: RingSecretKeyType,
        poly_length: usize,
        dimension: usize,
        modulus_minus_one: T,
        gaussian: Option<DiscreteGaussian<T>>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        let key_len = poly_length * dimension;
        let mut key = PolynomialOwned::zero(key_len);
        let distr = secret_key_type;
        match distr {
            RingSecretKeyType::Binary => key.random_binary_assign(rng),
            RingSecretKeyType::Ternary => key.random_ternary_assign(modulus_minus_one, rng),
            RingSecretKeyType::Gaussian => key.random_gaussian_assign(&gaussian.unwrap(), rng),
        };

        Self {
            key: key.into_owned(),
            poly_length,
            dimension,
            distr,
        }
    }
}

/// Represents a secret key for the (NTT) Ring Learning with Errors (RLWE) cryptographic scheme.
#[derive(Clone)]
pub struct NttGlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    poly_length: usize,
    dimension: usize,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> NttGlweSecretKey<T> {
    /// Creates a new [`NttGlweSecretKey<T>`].
    #[inline]
    pub fn new(
        key: Vec<T>,
        poly_length: usize,
        dimension: usize,
        distr: RingSecretKeyType,
    ) -> Self {
        Self {
            key,
            poly_length,
            dimension,
            distr,
        }
    }

    /// Returns the distr of this [`NttGlweSecretKey<T>`].
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    /// Creates a new [`NttGlweSecretKey`] from [`GlweSecretKey`].
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &GlweSecretKey<T>, ntt_table: &Table) -> Self
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = secret_key.poly_length;

        let mut key = secret_key.key.clone();
        key.chunks_exact_mut(poly_length)
            .for_each(|poly| ntt_table.transform_slice(poly));

        Self {
            key,
            poly_length,
            dimension: secret_key.dimension,
            distr: secret_key.distr,
        }
    }

    /// Performs `b-as`.
    pub fn phase_inplace<Table, M, S>(
        &self,
        cipher: &NttGlweCiphertext<S>,
        result: &mut PolynomialOwned<T>,
        ntt_table: &Table,
        modulus: M,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        S: RawData<Elem = T> + Data,
    {
        let mid = self.poly_length * self.dimension;
        let (a, b) = cipher.a_b_slices(mid);

        result.set_zero();

        a.chunks_exact(self.poly_length).for_each(|ai| {
            NttPolynomial(ArrayBase(result.as_mut())).add_mul_assign(
                &NttPolynomial(ArrayBase(ai)),
                &NttPolynomial(ArrayBase(self.key.as_ref())),
                modulus,
            );
        });
        NttPolynomial(ArrayBase(b))
            .sub_to_right(&mut NttPolynomial(ArrayBase(result.as_mut())), modulus);

        ntt_table.inverse_transform_slice(result.as_mut())
    }
}
