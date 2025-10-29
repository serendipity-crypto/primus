use primus_integer::UnsignedInteger;
use primus_lattice::rlwe::NttRlwe;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, Polynomial, RawData};
use primus_reduce::FieldContext;

use crate::{NttRlweCiphertext, NttRlweSecretKey, RlweParameters};

#[derive(Clone)]
pub struct NttRlwePublicKey<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    key: NttRlwe<S>,
}

impl<S, T> From<NttRlwe<S>> for NttRlwePublicKey<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    #[inline]
    fn from(key: NttRlwe<S>) -> Self {
        Self { key }
    }
}

impl<T: UnsignedInteger> NttRlwePublicKey<Vec<T>, T> {
    pub fn new<M, Table, R>(
        secret_key: &NttRlweSecretKey<T>,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) -> NttRlwePublicKey<Vec<T>>
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: rand::Rng + rand::CryptoRng,
    {
        let poly_length = params.poly_length();

        let mut data = vec![T::ZERO; poly_length * 2];

        let (a, b) = data.split_at_mut(poly_length);

        primus_distr::sample_gaussian_values_inplace(b, params.noise_distribution(), rng);

        ntt_table.transform_slice(b);

        primus_distr::sample_uniform_values_inplace(a, &params.cipher_modulus_uniform_distr(), rng);

        NttPolynomial(ArrayBase(b)).add_mul_assign(
            &NttPolynomial(ArrayBase(a)),
            secret_key,
            params.cipher_modulus(),
        );

        Self {
            key: NttRlwe::new(ArrayBase(data)),
        }
    }
}

impl<S, T> NttRlwePublicKey<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwePublicKey<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        Self {
            key: NttRlwe::from_bytes(data),
        }
    }
}

impl<S, T> NttRlwePublicKey<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwePublicKey<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        self.key.from_bytes_assign(data);
    }
}

impl<S, T> NttRlwePublicKey<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`NttRlwePublicKey<T>`] into bytes.
    #[inline]
    pub fn into_bytes(&self) -> Vec<u8> {
        self.key.to_bytes()
    }

    /// Converts [`NttRlwePublicKey<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn into_bytes_inplace(&self, data: &mut [u8]) {
        self.key.to_bytes_inplace(data);
    }

    /// Returns the bytes count of [`NttRlwePublicKey<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.key.bytes_count()
    }

    pub fn encrypt_inplace<M, Table, R, A, B>(
        &self,
        msg: &Polynomial<A>,
        result: &mut NttRlweCiphertext<B>,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: rand::Rng + rand::CryptoRng,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let modulus = params.cipher_modulus();

        let (a_in, b_in) = self.key.a_b_slices();
        let (a_out, b_out) = result.a_b_mut_slices();

        let mut v = vec![T::ZERO; poly_length];
        primus_distr::sample_ternary_values_inplace(&mut v, params.cipher_modulus_minus_one(), rng);
        ntt_table.transform_slice(&mut v);

        primus_distr::sample_gaussian_values_inplace(a_out, params.noise_distribution(), rng);
        ntt_table.transform_slice(a_out);
        NttPolynomial(ArrayBase(a_out)).add_mul_assign(
            &NttPolynomial(ArrayBase(a_in)),
            &NttPolynomial(ArrayBase(v.as_ref())),
            modulus,
        );

        primus_distr::sample_gaussian_values_inplace(b_out, params.noise_distribution(), rng);
        Polynomial(ArrayBase(&mut *b_out)).add_mul_scalar_assign(
            msg,
            params.plain_modulus_value(),
            modulus,
        );
        ntt_table.transform_slice(b_out);
        NttPolynomial(ArrayBase(b_out)).add_mul_assign(
            &NttPolynomial(ArrayBase(b_in)),
            &NttPolynomial(ArrayBase(v.as_ref())),
            modulus,
        );
    }

    pub fn encrypt<M, Table, R, A>(
        &self,
        msg: &Polynomial<A>,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) -> NttRlweCiphertext<Vec<T>>
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: rand::Rng + rand::CryptoRng,
        A: RawData<Elem = T> + Data,
    {
        let mut result = NttRlweCiphertext::zero(params.poly_length() * 2);
        self.encrypt_inplace(msg, &mut result, params, ntt_table, rng);
        result
    }

    pub fn encrypt_zeros_inplace<M, Table, R, A>(
        &self,
        result: &mut NttRlweCiphertext<A>,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: rand::Rng + rand::CryptoRng,
        A: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let modulus = params.cipher_modulus();

        let (a_in, b_in) = self.key.a_b_slices();
        let (a_out, b_out) = result.a_b_mut_slices();

        let mut v = vec![T::ZERO; poly_length];
        primus_distr::sample_ternary_values_inplace(&mut v, params.cipher_modulus_minus_one(), rng);
        ntt_table.transform_slice(&mut v);

        primus_distr::sample_gaussian_values_inplace(a_out, params.noise_distribution(), rng);
        ntt_table.transform_slice(a_out);
        NttPolynomial(ArrayBase(a_out)).add_mul_assign(
            &NttPolynomial(ArrayBase(a_in)),
            &NttPolynomial(ArrayBase(v.as_ref())),
            modulus,
        );

        primus_distr::sample_gaussian_values_inplace(b_out, params.noise_distribution(), rng);
        ntt_table.transform_slice(b_out);
        NttPolynomial(ArrayBase(b_out)).add_mul_assign(
            &NttPolynomial(ArrayBase(b_in)),
            &NttPolynomial(ArrayBase(v.as_ref())),
            modulus,
        );
    }

    pub fn encrypt_zeros<M, Table, R>(
        &self,
        params: &RlweParameters<T, M>,
        ntt_table: &Table,
        rng: &mut R,
    ) -> NttRlweCiphertext<Vec<T>>
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: rand::Rng + rand::CryptoRng,
    {
        let mut result = NttRlweCiphertext::zero(params.poly_length() * 2);
        self.encrypt_zeros_inplace(&mut result, params, ntt_table, rng);
        result
    }
}
