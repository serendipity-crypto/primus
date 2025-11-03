use primus_integer::{UnsignedInteger, izip};
use primus_lattice::{ggsw::DcrtGgsw, glwe::DcrtGlwe};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use crate::{CrtGgswParameters, CrtGlevParameters, CrtGlweParameters, DcrtGlweSecretKey};

#[derive(Clone)]
pub struct DcrtGlwePublicKey<T: UnsignedInteger> {
    key: DcrtGlwe<Vec<T>>,
}

impl<T: UnsignedInteger> DcrtGlwePublicKey<T> {
    /// Creates a new [`DcrtGlwePublicKey<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        self.key.from_bytes_assign(data);
    }

    /// Converts [`DcrtGlwePublicKey<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        self.key.to_bytes_inplace(data);
    }

    /// Returns the bytes count of [`DcrtGlwePublicKey<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.key.bytes_count()
    }

    pub fn new<Table, R, M>(
        secret_key: &DcrtGlweSecretKey<T>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) -> DcrtGlwePublicKey<T>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let rns_glwe_len = params.rns_glwe_len();

        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let mut data: DcrtGlwe<Vec<T>> = DcrtGlwe::zero(rns_glwe_len);

        let (a, b) = data.a_b_mut_slices(params.rns_glwe_mid());

        primus_distr::sample_crt_gaussian_values_inplace(
            b,
            poly_length,
            moduli_value,
            params.noise_distribution(),
            rng,
        );

        table.transform_slice(b);
        let mut b_dcrt_poly = DcrtPolynomial(ArrayBase(b));

        a.chunks_exact_mut(rns_poly_len)
            .zip(secret_key.iter_dcrt_poly())
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

        Self { key: data }
    }

    pub fn encrypt<R, M, Table, A>(
        &self,
        message: &CrtPolynomial<A>,
        params: &CrtGlweParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) -> DcrtGlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
    {
        let dimension = params.dimension();
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let rns_glwe_len = params.rns_glwe_len();

        let moduli = params.cipher_moduli();

        let mut result = DcrtGlwe::new(ArrayBase(vec![T::ZERO; rns_glwe_len]));
        let mut v = vec![T::ZERO; rns_poly_len];

        primus_distr::sample_crt_ternary_values_inplace(
            &mut v,
            poly_length,
            params.cipher_moduli_minus_one(),
            rng,
        );
        table.transform_slice(&mut v);
        let v_dcrt_poly = DcrtPolynomial(ArrayBase(v.as_ref()));

        result
            .iter_dcrt_poly_mut(rns_poly_len)
            .zip(self.key.iter_dcrt_poly(rns_poly_len))
            .enumerate()
            .for_each(|(i, (ai, pk_ai))| {
                primus_distr::sample_crt_gaussian_values_inplace(
                    ai,
                    poly_length,
                    params.cipher_moduli_value(),
                    params.noise_distribution(),
                    rng,
                );

                if i == dimension {
                    CrtPolynomial(ArrayBase(&mut *ai)).add_assign(message, poly_length, moduli);
                }
                table.transform_slice(ai);
                DcrtPolynomial(ArrayBase(ai)).add_mul_assign(
                    &v_dcrt_poly,
                    &DcrtPolynomial(ArrayBase(pk_ai)),
                    poly_length,
                    moduli,
                );
            });

        result
    }

    fn encrypt_monomial_in_dcrt_glev_inplace<R, Table, M>(
        &self,
        index: usize,
        coeff_residues: &[T],
        degree: usize,
        dcrt_glev: &mut [T],
        params: &CrtGlevParameters<T, M>,
        table: &Table,
        v: &mut [T],
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let rns_glwe_len = params.rns_glwe_len();

        let moduli = params.cipher_moduli();

        izip!(
            dcrt_glev.chunks_exact_mut(rns_glwe_len),
            params.basis().iter_scalar_residues()
        )
        .for_each(|(dcrt_glwe, scalar_residues)| {
            primus_distr::sample_crt_ternary_values_inplace(
                v,
                poly_length,
                params.cipher_moduli_minus_one(),
                rng,
            );
            table.transform_slice(v);
            let v_poly = DcrtPolynomial(ArrayBase(v.as_ref()));

            dcrt_glwe
                .chunks_exact_mut(rns_poly_len)
                .zip(self.key.iter_dcrt_poly(rns_poly_len))
                .enumerate()
                .for_each(|(i, (ai, pk_ai))| {
                    primus_distr::sample_crt_gaussian_values_inplace(
                        ai,
                        poly_length,
                        params.cipher_moduli_value(),
                        params.noise_distribution(),
                        rng,
                    );

                    if i == index {
                        izip!(
                            ai.chunks_exact_mut(poly_length),
                            coeff_residues,
                            scalar_residues,
                            moduli,
                        )
                        .for_each(
                            |(poly, &coeff_residue, &scalar_residue, &modulus)| {
                                poly[degree] = modulus.reduce_mul_add(
                                    coeff_residue,
                                    scalar_residue,
                                    poly[degree],
                                );
                            },
                        );
                    }

                    table.transform_slice(ai);

                    DcrtPolynomial(ArrayBase(ai)).add_mul_assign(
                        &v_poly,
                        &DcrtPolynomial(ArrayBase(pk_ai)),
                        poly_length,
                        moduli,
                    );
                });
        });
    }

    /// Generate a [`DcrtGgsw`] ciphertext which encrypted `coeff*X^degree`.
    pub fn encrypt_monomial_ggsw<R, Table, M>(
        &self,
        coeff_residues: &[T],
        degree: usize,
        params: &CrtGgswParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) -> DcrtGgsw<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        M: FieldContext<T>,
    {
        let rns_ggsw_len = params.rns_ggsw_len();

        let mut result: DcrtGgsw<Vec<T>> = DcrtGgsw::zero(rns_ggsw_len);

        self.encrypt_monomial_ggsw_inplace(coeff_residues, degree, &mut result, params, table, rng);

        result
    }

    /// Generate a [`DcrtGgsw`] ciphertext which encrypted `coeff*X^degree`.
    pub fn encrypt_monomial_ggsw_inplace<R, Table, M, A>(
        &self,
        coeff_residues: &[T],
        degree: usize,
        result: &mut DcrtGgsw<A>,
        params: &CrtGgswParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        let rns_poly_len = params.rns_poly_len();
        let rns_glev_len = params.rns_glev_len();

        assert_eq!(result.as_ref().len(), params.rns_ggsw_len());
        assert!(degree < params.poly_length());

        let mut v = vec![T::ZERO; rns_poly_len];

        result
            .iter_dcrt_glev_mut(rns_glev_len)
            .enumerate()
            .for_each(|(i, dcrt_glev)| {
                self.encrypt_monomial_in_dcrt_glev_inplace(
                    i,
                    coeff_residues,
                    degree,
                    dcrt_glev,
                    params,
                    table,
                    &mut v,
                    rng,
                );
            });
    }
}
