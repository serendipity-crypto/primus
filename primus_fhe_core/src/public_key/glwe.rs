use primus_integer::{UnsignedInteger, izip};
use primus_lattice::{ggsw::DcrtGgsw, glwe::DcrtGlwe};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use crate::{CrtGgswParameters, CrtGlevParameters, CrtGlweParameters, DcrtGlweSecretKey};

#[derive(Clone)]
pub struct DcrtGlwePublicKey<T: UnsignedInteger> {
    key: DcrtGlwe<Vec<T>>,
    moduli_count: usize,
    poly_length: usize,
    dimension: usize,
    crt_poly_length: usize, // moduli_count * poly_length
    crt_glwe_len: usize,    // moduli_count * poly_length * (dimension + 1)
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

    pub fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    pub fn dimension(&self) -> usize {
        self.dimension
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
        let poly_length = secret_key.poly_length();
        let dimension = secret_key.dimension();
        let crt_poly_length = secret_key.crt_poly_length();
        let crt_glwe_len = secret_key.crt_glwe_len();

        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let mut data: DcrtGlwe<Vec<T>> = DcrtGlwe::zero(crt_glwe_len);

        let (a, b) = data.a_b_mut_slices(crt_glwe_len - crt_poly_length);

        a.chunks_exact_mut(crt_poly_length).for_each(|ai| {
            primus_distr::sample_crt_uniform_values_inplace(ai, poly_length, uniform_distrs, rng);
        });
        primus_distr::sample_crt_gaussian_values_inplace(
            b,
            poly_length,
            moduli_value,
            params.noise_distribution(),
            rng,
        );

        table.transform_slice(b);

        let mut b_poly = DcrtPolynomial(ArrayBase(b));

        a.chunks_exact_mut(crt_poly_length)
            .zip(secret_key.iter_dcrt_poly())
            .for_each(|(ai, si)| {
                b_poly.add_mul_assign(
                    &DcrtPolynomial(ArrayBase(ai)),
                    &DcrtPolynomial(ArrayBase(si)),
                    poly_length,
                    moduli,
                );
            });

        Self {
            key: data,
            moduli_count: secret_key.moduli_count(),
            poly_length,
            dimension,
            crt_poly_length,
            crt_glwe_len,
        }
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
        let poly_length = self.poly_length;
        let crt_glwe_len = self.crt_glwe_len;
        let moduli = params.cipher_moduli();

        let mut result = DcrtGlwe::new(ArrayBase(vec![T::ZERO; crt_glwe_len]));
        let mut v = vec![T::ZERO; self.crt_poly_length];

        primus_distr::sample_crt_ternary_values_inplace(
            &mut v,
            poly_length,
            params.cipher_moduli_minus_one(),
            rng,
        );
        table.transform_slice(&mut v);
        let v_poly = DcrtPolynomial(ArrayBase(v.as_ref()));

        result
            .iter_dcrt_poly_mut(self.crt_poly_length)
            .zip(self.key.iter_dcrt_poly(self.crt_poly_length))
            .enumerate()
            .for_each(|(i, (ai, pk_ai))| {
                primus_distr::sample_crt_gaussian_values_inplace(
                    ai,
                    poly_length,
                    params.cipher_moduli_value(),
                    params.noise_distribution(),
                    rng,
                );

                if i == self.dimension {
                    CrtPolynomial(ArrayBase(&mut *ai)).add_assign(message, poly_length, moduli);
                }
                table.transform_slice(ai);
                DcrtPolynomial(ArrayBase(ai)).add_mul_assign(
                    &v_poly,
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
        let poly_length = self.poly_length;
        let crt_poly_length = self.crt_poly_length;
        let crt_glwe_len = self.crt_glwe_len;
        let moduli = params.cipher_moduli();

        izip!(
            dcrt_glev.chunks_exact_mut(crt_glwe_len),
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
                .chunks_exact_mut(crt_poly_length)
                .zip(self.key.iter_dcrt_poly(crt_poly_length))
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
        let dcrt_glwe_len = self.crt_glwe_len;
        let dcrt_glev_len = params.basis().decompose_length() * dcrt_glwe_len;
        let dcrt_ggsw_len = (self.dimension + 1) * dcrt_glev_len;

        let mut dcrt_ggsw: DcrtGgsw<Vec<T>> = DcrtGgsw::zero(dcrt_ggsw_len);

        let mut v = vec![T::ZERO; self.crt_poly_length];

        dcrt_ggsw
            .iter_dcrt_glev_mut(dcrt_glev_len)
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

        dcrt_ggsw
    }
}
