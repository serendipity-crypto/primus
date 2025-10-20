use primus_integer::UnsignedInteger;
use primus_ntt::{Dcrt, DcrtTable, Ntt, NttTable};
use primus_poly::{
    ArrayBase, Data, DataMut, NttPolynomial, PolynomialOwned, RawData, crt::CrtPolynomial,
    dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::{CrtGlweParameters, DcrtGlweCiphertext, GlweParameters, NttGlweCiphertext};

use super::RingSecretKeyType;

/// Represents a secret key for the Module Learning with Errors (MLWE) cryptographic scheme.
#[derive(Clone)]
pub struct GlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    poly_length: usize,
    dimension: usize,
    distr: RingSecretKeyType,
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
        debug_assert!(poly_length.is_power_of_two());
        debug_assert_eq!(key.len(), poly_length * dimension);
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
    pub fn generate<R, M>(params: &GlweParameters<T, M>, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let dimension = params.dimension();

        let key_len = poly_length * dimension;
        let mut key = PolynomialOwned::zero(key_len);
        let distr = params.secret_key_type();
        match distr {
            RingSecretKeyType::Binary => key.random_binary_assign(rng),
            RingSecretKeyType::Ternary => {
                key.random_ternary_assign(params.cipher_modulus_minus_one(), rng)
            }
            RingSecretKeyType::Gaussian => {
                // FIXME
                key.random_gaussian_assign(params.noise_distribution(), rng)
            }
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
        debug_assert!(poly_length.is_power_of_two());
        debug_assert_eq!(key.len(), poly_length * dimension);
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

        let mut result_poly = NttPolynomial(ArrayBase(result.as_mut()));

        a.chunks_exact(self.poly_length).for_each(|ai| {
            result_poly.add_mul_assign(
                &NttPolynomial(ArrayBase(ai)),
                &NttPolynomial(ArrayBase(self.key.as_ref())),
                modulus,
            );
        });
        NttPolynomial(ArrayBase(b)).sub_to_right(&mut result_poly, modulus);

        ntt_table.inverse_transform_slice(result.as_mut())
    }
}

pub struct CrtGlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    moduli_count: usize,
    poly_length: usize,
    dimension: usize,
    crt_poly_length: usize,
    distr: RingSecretKeyType,
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
        debug_assert!(poly_length.is_power_of_two());
        let crt_poly_length = moduli_count * poly_length;
        debug_assert_eq!(key.len(), crt_poly_length * dimension);
        Self {
            key,
            moduli_count,
            poly_length,
            dimension,
            crt_poly_length,
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

    /// Returns the crt poly length of this [`CrtGlweSecretKey<T>`].
    pub fn crt_poly_length(&self) -> usize {
        self.crt_poly_length
    }

    pub fn iter_crt_poly(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.crt_poly_length)
    }

    /// Returns the distr of this [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    pub fn generate<R, M>(params: &CrtGlweParameters<T, M>, rng: &mut R) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let dimension = params.dimension();
        let moduli_value = params.cipher_moduli_value();
        let moduli_minus_one = params.cipher_moduli_minus_one();
        let moduli_count = params.cipher_moduli_count();

        let crt_poly_length = moduli_count * poly_length;

        let distr = params.secret_key_type();

        let mut key = vec![T::ZERO; crt_poly_length * dimension];

        match distr {
            RingSecretKeyType::Binary => {
                key.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
                    primus_distr::sample_crt_binary_values_inplace(crt_poly, poly_length, rng);
                });
            }
            RingSecretKeyType::Ternary => {
                key.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
                    primus_distr::sample_crt_ternary_values_inplace(
                        crt_poly,
                        poly_length,
                        moduli_minus_one,
                        rng,
                    );
                });
            }
            RingSecretKeyType::Gaussian => {
                // FIXME
                key.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
                    primus_distr::sample_crt_gaussian_values_inplace(
                        crt_poly,
                        poly_length,
                        moduli_value,
                        params.noise_distribution(),
                        rng,
                    );
                });
            }
        };

        Self {
            key,
            moduli_count,
            poly_length,
            dimension,
            crt_poly_length,
            distr,
        }
    }
}

pub struct DcrtGlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    moduli_count: usize,
    poly_length: usize,
    dimension: usize,
    crt_poly_length: usize,
    crt_glwe_len: usize,
    distr: RingSecretKeyType,
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

    /// Returns the crt poly length of this [`DcrtGlweSecretKey<T>`].
    pub fn crt_poly_length(&self) -> usize {
        self.crt_poly_length
    }

    /// Returns the crt glwe len of this [`DcrtGlweSecretKey<T>`].
    pub fn crt_glwe_len(&self) -> usize {
        self.crt_glwe_len
    }

    /// Returns the distr of this [`DcrtGlweSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    pub fn iter_dcrt_poly(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.crt_poly_length)
    }

    /// Creates a new [`NttGlweSecretKey`] from [`GlweSecretKey`].
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &CrtGlweSecretKey<T>, table: &Table) -> Self
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let moduli_count = secret_key.moduli_count;
        let poly_length = secret_key.poly_length;
        let dimension = secret_key.dimension;
        let crt_poly_length = secret_key.crt_poly_length;

        let mut key = secret_key.key.clone();

        key.chunks_exact_mut(crt_poly_length).for_each(|crt_poly| {
            table.transform_slice(crt_poly);
        });

        Self {
            key,
            moduli_count,
            poly_length,
            dimension,
            crt_poly_length,
            crt_glwe_len: crt_poly_length * (dimension + 1),
            distr: secret_key.distr,
        }
    }

    pub fn encrypt_inplace<R, M, Table, A, B>(
        &self,
        msg: &CrtPolynomial<A>,
        result: &mut DcrtGlweCiphertext<B>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let crt_poly_length = self.crt_poly_length;
        let moduli = params.cipher_moduli();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let (a, b) = result.a_b_mut_slices(self.crt_glwe_len - crt_poly_length);

        primus_distr::sample_crt_gaussian_values_inplace(
            b,
            poly_length,
            params.cipher_moduli_value(),
            params.noise_distribution(),
            rng,
        );

        let mut b_crt_poly = CrtPolynomial(ArrayBase(b));
        b_crt_poly.add_mul_scalar_residues_assign(
            msg,
            params.delta_residues(),
            poly_length,
            moduli,
        );
        let mut b_dcrt_poly = table.transform_inplace(b_crt_poly);

        a.chunks_exact_mut(crt_poly_length)
            .zip(self.iter_dcrt_poly())
            .for_each(|(ai, si)| {
                primus_distr::sample_crt_uniform_values_inplace(
                    ai,
                    poly_length,
                    uniform_distrs,
                    rng,
                );

                b_dcrt_poly.add_mul_assign(
                    &DcrtPolynomial(ArrayBase(ai)),
                    &DcrtPolynomial(ArrayBase(si)),
                    poly_length,
                    moduli,
                );
            });
    }

    /// Performs `b-as`.
    pub fn phase_inplace<M, Table, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        msg: &mut CrtPolynomial<B>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let moduli = params.cipher_moduli();

        let mut temp = DcrtPolynomial::new(ArrayBase(msg.as_mut()));
        temp.set_zero();

        let (a, b) = ciphertext.a_b_slices(self.crt_glwe_len - self.crt_poly_length);

        // temp = ∑a*s
        a.chunks_exact(self.crt_poly_length)
            .zip(self.iter_dcrt_poly())
            .for_each(|(ai, si)| {
                temp.add_mul_assign(
                    &DcrtPolynomial(ArrayBase(ai)),
                    &DcrtPolynomial(ArrayBase(si)),
                    poly_length,
                    moduli,
                );
            });

        // temp = b - ∑a*s
        DcrtPolynomial(ArrayBase(b)).sub_to_right(&mut temp, poly_length, moduli);
        table.inverse_transform_slice(temp.as_mut());
    }

    pub fn decrypt_inplace<M, Table, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        msg: &mut CrtPolynomial<B>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let inverse_delta_residues = params.inverse_delta_residues();
        self.phase_inplace(ciphertext, msg, params, table);
        msg.mul_scalar_residues_assign(
            inverse_delta_residues,
            self.poly_length,
            params.cipher_moduli(),
        );
    }
}
