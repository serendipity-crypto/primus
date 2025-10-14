use primus_distr::{DiscreteGaussian, SignedDiscreteGaussian};
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt, NttTable};
use primus_poly::{ArrayBase, Data, NttPolynomial, PolynomialOwned, RawData, crt::CrtPolynomial};
use primus_reduce::FieldContext;

use crate::{DcrtGlweCiphertext, NttGlweCiphertext};

use super::RingSecretKeyType;

/// Represents a secret key for the Module Learning with Errors (MLWE) cryptographic scheme.
#[derive(Clone)]
pub struct GlweSecretKey<T: UnsignedInteger> {
    pub(crate) key: Vec<T>,
    pub(crate) poly_length: usize,
    pub(crate) dimension: usize,
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

    /// Returns the dimension of this [`GlweSecretKey<T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
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

    /// Returns the poly length of this [`NttGlweSecretKey<T>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the dimension of this [`NttGlweSecretKey<T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
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

pub struct CrtGlweSecretKey<T: UnsignedInteger> {
    pub(crate) key: Vec<T>,
    pub(crate) moduli_count: usize,
    pub(crate) poly_length: usize,
    pub(crate) dimension: usize,
    pub(crate) distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> CrtGlweSecretKey<T> {
    /// Creates a new [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn new(
        key: Vec<T>,
        moduli_count: usize,
        poly_length: usize,
        dimension: usize,
        distr: RingSecretKeyType,
    ) -> Self {
        Self {
            key,
            moduli_count,
            poly_length,
            dimension,
            distr,
        }
    }

    /// Returns the moduli count of this [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    /// Returns the poly length of this [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the dimension of this [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the distr of this [`CrtGlweSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    pub fn generate<R>(
        secret_key_type: RingSecretKeyType,
        poly_length: usize,
        dimension: usize,
        moduli_minus_one: &[T],
        gaussian: Option<&SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        let moduli_count = moduli_minus_one.len();

        let key = match secret_key_type {
            RingSecretKeyType::Binary => {
                primus_distr::sample_crt_binary_values(poly_length * dimension, moduli_count, rng)
            }
            RingSecretKeyType::Ternary => primus_distr::sample_crt_ternary_values(
                poly_length * dimension,
                moduli_minus_one,
                rng,
            ),
            RingSecretKeyType::Gaussian => {
                let moduli: Vec<T> = moduli_minus_one.iter().map(|&v| v + T::ONE).collect();
                primus_distr::sample_crt_gaussian_values(
                    poly_length * dimension,
                    &moduli,
                    gaussian.unwrap(),
                    rng,
                )
            }
        };

        Self {
            key,
            moduli_count,
            poly_length,
            dimension,
            distr: secret_key_type,
        }
    }
}

pub struct DcrtGlweSecretKey<T: UnsignedInteger> {
    pub(crate) key: Vec<T>,
    pub(crate) moduli_count: usize,
    pub(crate) poly_length: usize,
    pub(crate) dimension: usize,
    pub(crate) distr: RingSecretKeyType,
    single_modulus_len: usize,
}

impl<T: UnsignedInteger> DcrtGlweSecretKey<T> {
    /// Returns the moduli count of this [`DcrtGlweSecretKey<T>`].
    pub fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    /// Returns the poly length of this [`DcrtGlweSecretKey<T>`].
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the dimension of this [`DcrtGlweSecretKey<T>`].
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the distr of this [`DcrtGlweSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    #[inline]
    pub fn iter_each_modulus(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.single_modulus_len)
    }

    /// Creates a new [`NttGlweSecretKey`] from [`GlweSecretKey`].
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &CrtGlweSecretKey<T>, table: &Table) -> Self
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = secret_key.poly_length;
        let dimension = secret_key.dimension;
        let single_modulus_len = poly_length * dimension;

        let mut key = secret_key.key.clone();

        key.chunks_exact_mut(single_modulus_len)
            .zip(table.iter())
            .for_each(|(chunk, ntt_table)| {
                chunk.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
            });

        Self {
            key,
            moduli_count: secret_key.moduli_count,
            poly_length,
            dimension,
            distr: secret_key.distr,
            single_modulus_len,
        }
    }

    /// Performs `b-as`.
    pub fn phase_inplace<Table, M, S>(
        &self,
        cipher: &DcrtGlweCiphertext<S>,
        result: &mut CrtPolynomial<Vec<T>>,
        table: &Table,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        S: RawData<Elem = T> + Data,
    {
        let mid = self.poly_length * self.dimension;
        let single_modulus_len = self.poly_length + mid;

        result.set_zero();

        izip!(
            cipher.data.chunks_exact(single_modulus_len),
            result.iter_each_modulus_mut(self.poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, poly, ntt_table, modulus)| {
            let (a, b) = glwe.split_at(mid);
            let mut res = NttPolynomial(ArrayBase(&mut *poly));

            izip!(
                a.chunks_exact(self.poly_length),
                self.key.chunks_exact(self.poly_length)
            )
            .for_each(|(a, s)| {
                res.add_mul_assign(
                    &NttPolynomial(ArrayBase(a)),
                    &NttPolynomial(ArrayBase(s)),
                    *modulus,
                );
            });

            NttPolynomial(ArrayBase(b)).sub_to_right(&mut res, *modulus);

            ntt_table.inverse_transform_slice(poly);
        });
    }
}
