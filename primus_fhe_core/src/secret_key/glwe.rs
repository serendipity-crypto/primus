use itertools::Itertools;
use primus_factor::FactorMul;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::glev::DcrtGlev;
use primus_ntt::{Dcrt, DcrtTable, Ntt, NttTable};
use primus_poly::{
    ArrayBase, Data, DataMut, NttPolynomial, Polynomial, PolynomialOwned, RawData,
    crt::CrtPolynomial, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::{
    CrtGlevParameters, CrtGlweParameters, DcrtGlweCiphertext, GlweParameters, NttGlweCiphertext,
};

use super::RingSecretKeyType;

/// Represents a secret key for the Module Learning with Errors (MLWE) cryptographic scheme.
#[derive(Clone)]
pub struct GlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    dimension: usize,
    poly_length: usize,
    distr: RingSecretKeyType,
}

impl<T: UnsignedInteger> GlweSecretKey<T> {
    /// Creates a new [`GlweSecretKey<T>`].
    #[inline]
    pub fn new(
        key: Vec<T>,
        dimension: usize,
        poly_length: usize,
        distr: RingSecretKeyType,
    ) -> Self {
        debug_assert!(poly_length.is_power_of_two());
        debug_assert_eq!(key.len(), poly_length * dimension);
        Self {
            key,
            dimension,
            poly_length,
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
        let dimension = params.dimension();
        let poly_length = params.poly_length();

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

#[derive(Clone)]
pub struct CrtGlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    distr: RingSecretKeyType,
    rns_poly_len: usize,
}

impl<T: UnsignedInteger> CrtGlweSecretKey<T> {
    /// Creates a new [`CrtGlweSecretKey<T>`].
    #[inline]
    pub fn new(key: Vec<T>, distr: RingSecretKeyType, rns_poly_len: usize) -> Self {
        Self {
            key,
            distr,
            rns_poly_len,
        }
    }

    pub fn key(&self) -> &[T] {
        &self.key
    }

    pub fn key_mut(&mut self) -> &mut Vec<T> {
        &mut self.key
    }

    pub fn iter_crt_poly(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.rns_poly_len)
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
        let rns_poly_len = params.rns_poly_len();

        let secret_key_type = params.secret_key_type();

        let mut key = vec![T::ZERO; params.secret_key_len()];

        match secret_key_type {
            RingSecretKeyType::Binary => {
                key.chunks_exact_mut(rns_poly_len).for_each(|crt_poly| {
                    primus_distr::sample_crt_binary_values_inplace(crt_poly, poly_length, rng);
                });
            }
            RingSecretKeyType::Ternary => {
                let moduli_minus_one = params.cipher_moduli_minus_one();
                key.chunks_exact_mut(rns_poly_len).for_each(|crt_poly| {
                    primus_distr::sample_crt_ternary_values_inplace(
                        crt_poly,
                        poly_length,
                        moduli_minus_one,
                        rng,
                    );
                });
            }
            RingSecretKeyType::Gaussian => {
                unimplemented!()
            }
        };

        Self {
            key,
            distr: secret_key_type,
            rns_poly_len,
        }
    }
}

pub struct DcrtGlweSecretKey<T: UnsignedInteger> {
    key: Vec<T>,
    distr: RingSecretKeyType,
    rns_poly_len: usize,
}

impl<T: UnsignedInteger> DcrtGlweSecretKey<T> {
    pub fn key(&self) -> &[T] {
        &self.key
    }

    /// Returns the distr of this [`DcrtGlweSecretKey<T>`].
    pub fn distr(&self) -> RingSecretKeyType {
        self.distr
    }

    pub fn iter_dcrt_poly(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.rns_poly_len)
    }

    /// Creates a new [`NttGlweSecretKey`] from [`GlweSecretKey`].
    #[inline]
    pub fn from_coeff_secret_key<Table>(secret_key: &CrtGlweSecretKey<T>, table: &Table) -> Self
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let rns_poly_len = secret_key.rns_poly_len;

        let mut key = secret_key.key.clone();

        key.chunks_exact_mut(rns_poly_len).for_each(|crt_poly| {
            table.transform_slice(crt_poly);
        });

        Self {
            key,
            distr: secret_key.distr,
            rns_poly_len,
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
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let rns_glwe_mid = params.rns_glwe_mid();
        let moduli = params.cipher_moduli();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let (a, b) = result.a_b_mut_slices(rns_glwe_mid);

        primus_distr::sample_crt_gaussian_values_inplace(
            b,
            poly_length,
            params.cipher_moduli_value(),
            params.noise_distribution(),
            rng,
        );

        let mut b_crt_poly = CrtPolynomial(ArrayBase(b));
        b_crt_poly.add_mul_scalar_assign(msg, params.delta_mod_q(), poly_length, moduli);
        let mut b_dcrt_poly = table.transform_inplace(b_crt_poly);

        a.chunks_exact_mut(rns_poly_len)
            .zip_eq(self.iter_dcrt_poly())
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

    fn encrypt_custom_delta_dcrt_glwe_inplace<R, M, Table, A, B>(
        &self,
        msg: &CrtPolynomial<A>,
        delta_residues: &[T],
        result: &mut DcrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let moduli = params.cipher_moduli();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let (a, b) = result.a_b_mut_slices(params.rns_glwe_mid());

        primus_distr::sample_crt_gaussian_values_inplace(
            b,
            poly_length,
            params.cipher_moduli_value(),
            params.noise_distribution(),
            &mut *rng,
        );

        let mut b_crt_poly = CrtPolynomial(ArrayBase(b));
        b_crt_poly.add_mul_scalar_assign(&msg, delta_residues, poly_length, moduli);
        let mut b_dcrt_poly = table.transform_inplace(b_crt_poly);

        a.chunks_exact_mut(rns_poly_len)
            .zip(self.iter_dcrt_poly())
            .for_each(|(ai, si)| {
                primus_distr::sample_crt_uniform_values_inplace(
                    ai,
                    poly_length,
                    uniform_distrs,
                    &mut *rng,
                );
                b_dcrt_poly.add_mul_assign(
                    &DcrtPolynomial(ArrayBase(ai)),
                    &DcrtPolynomial(ArrayBase(si)),
                    poly_length,
                    moduli,
                );
            });
    }

    pub fn encrypt_dcrt_glev_inplace<R, M, Table, A, B>(
        &self,
        msg: &CrtPolynomial<A>,
        result: &mut DcrtGlev<B>,
        params: &CrtGlevParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        result
            .iter_dcrt_glwe_mut(params.rns_glwe_len())
            .zip(params.basis().iter_scalar_residues())
            .for_each(|(dcrt_glwe, scalar_residues)| {
                self.encrypt_custom_delta_dcrt_glwe_inplace(
                    msg,
                    scalar_residues,
                    &mut DcrtGlweCiphertext::new(ArrayBase(dcrt_glwe)),
                    params,
                    table,
                    rng,
                );
            });
    }

    /// Performs `b - ∑ a*s`.
    fn phase_inplace<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        msg_mod_q: &mut DcrtPolynomial<B>,
        params: &CrtGlweParameters<T, M>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let moduli = params.cipher_moduli();

        let (a, b) = ciphertext.a_b_slices(params.rns_glwe_mid());

        msg_mod_q.set_zero();

        // msg_mod_q = ∑a*s
        a.chunks_exact(params.rns_poly_len())
            .zip(self.iter_dcrt_poly())
            .for_each(|(ai, si)| {
                msg_mod_q.add_mul_assign(
                    &DcrtPolynomial(ArrayBase(ai)),
                    &DcrtPolynomial(ArrayBase(si)),
                    poly_length,
                    moduli,
                );
            });

        // msg_mod_q = b - ∑ a*s
        DcrtPolynomial(ArrayBase(b)).sub_to_right(msg_mod_q, poly_length, moduli);
    }

    pub fn decrypt<M, Table, A>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
        context: &mut DcrtGlweDecryptContext<T>,
    ) -> PolynomialOwned<T>
    where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
    {
        let mut msg = PolynomialOwned::zero(params.poly_length());
        self.decrypt_inplace(ciphertext, &mut msg, params, table, context);
        msg
    }

    pub fn decrypt_inplace<M, Table, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        msg: &mut Polynomial<B>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
        context: &mut DcrtGlweDecryptContext<T>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let q = params.cipher_moduli();
        let t = params.plain_modulus_value();
        let gamma = params.gamma();
        let inv_gamma_mod_t = params.inv_gamma_mod_t();

        let (msg_mod_q, msg_mod_t_gamma) = context.as_mut();

        self.phase_inplace(ciphertext, msg_mod_q, params);

        table.inverse_transform_slice(msg_mod_q.as_mut());

        msg_mod_q.mul_scalar_assign(params.t_gamma_mod_q(), poly_length, q);

        params.converter().fast_convert_array(
            msg_mod_q.as_ref(),
            msg_mod_t_gamma.as_mut(),
            poly_length,
        );

        msg_mod_t_gamma.mul_scalar_assign(
            params.minus_inv_q_mod_t_gamma(),
            poly_length,
            params.t_gamma(),
        );

        let (y_t_slices, y_gamma_slices) = msg_mod_t_gamma.as_ref().split_at(poly_length);

        izip!(msg.iter_mut(), y_t_slices, y_gamma_slices).for_each(|(res, &y_t, &y_gamma)| {
            let mut temp = gamma - y_gamma + y_t;
            if temp >= gamma {
                temp -= gamma;
            }
            *res = inv_gamma_mod_t.factor_mul_modulo(temp, t);
        });
    }
}

pub struct DcrtGlweDecryptContext<T: UnsignedInteger> {
    msg_mod_q: DcrtPolynomial<Vec<T>>,
    msg_mod_t_gamma: CrtPolynomial<Vec<T>>,
}

impl<T: UnsignedInteger> DcrtGlweDecryptContext<T> {
    /// Creates a new [`DcrtGlweDecryptContext<T>`].
    #[inline]
    pub fn new(moduli_count: usize, poly_length: usize) -> Self {
        let msg_mod_q: DcrtPolynomial<Vec<T>> = DcrtPolynomial::zero(moduli_count * poly_length);
        let msg_mod_t_gamma: CrtPolynomial<Vec<T>> = CrtPolynomial::zero(2 * poly_length);

        Self {
            msg_mod_q,
            msg_mod_t_gamma,
        }
    }

    #[inline]
    pub fn as_mut(
        &mut self,
    ) -> (
        &mut DcrtPolynomial<Vec<T>, T>,
        &mut CrtPolynomial<Vec<T>, T>,
    ) {
        (&mut self.msg_mod_q, &mut self.msg_mod_t_gamma)
    }
}
