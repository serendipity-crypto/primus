use std::slice::IterMut;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::{DcrtGgsw, DcrtGlwe};
use primus_modulo::AddModulo;
use primus_modulo::MulAddModulo;
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, NttPolynomial, Polynomial, RawData, crt::CrtPolynomial,
    dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use rand::distr::Uniform;

use crate::DcrtGlweSecretKey;

#[derive(Clone)]
pub struct DcrtGlwePublicKey<T: UnsignedInteger, M: FieldContext<T>> {
    key: DcrtGlwe<Vec<T>>,
    moduli_count: usize,
    poly_length: usize,
    dimension: usize,
    a_b_mid: usize,      // poly_length * dimension
    glwe_len: usize,     // poly_length * (dimension + 1)
    crt_glwe_len: usize, // moduli_count * poly_length * (dimension + 1)
    moduli: Vec<M>,
    modulus_values: Vec<T>,
}

impl<T: UnsignedInteger, M: FieldContext<T>> DcrtGlwePublicKey<T, M> {
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

    pub fn moduli(&self) -> &[M] {
        &self.moduli
    }

    pub fn modulus_values(&self) -> &[T] {
        &self.modulus_values
    }

    pub fn iter_each_modulus(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.as_ref().chunks_exact(self.glwe_len)
    }

    pub fn new<Table, R>(
        secret_key: &DcrtGlweSecretKey<T>,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        moduli: &[M],
        rng: &mut R,
    ) -> DcrtGlwePublicKey<T, M>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let moduli_count = secret_key.moduli_count;
        let poly_length = secret_key.poly_length;
        let dimension = secret_key.dimension;

        let a_b_mid = poly_length * dimension;
        let glwe_len = a_b_mid + poly_length;
        let key_len = moduli_count * glwe_len;

        let modulus_values: Vec<T> = moduli.iter().map(|m| m.value_unchecked()).collect();
        let uniform_distrs: Vec<Uniform<T>> =
            moduli.iter().map(|m| m.uniform_distribution()).collect();

        let mut data = vec![T::ZERO; key_len];

        let (a_iters, b_iters): (Vec<IterMut<'_, T>>, Vec<IterMut<'_, T>>) = data
            .chunks_exact_mut(glwe_len)
            .map(|glwe| {
                let (a, b) = unsafe { glwe.split_at_mut_unchecked(a_b_mid) };
                (a.iter_mut(), b.iter_mut())
            })
            .collect();

        primus_distr::sample_crt_uniform_values_iter_mut(a_iters, &uniform_distrs, rng);
        primus_distr::sample_crt_gaussian_values_iter_mut(b_iters, &modulus_values, gaussian, rng);

        izip!(
            data.chunks_exact_mut(glwe_len),
            secret_key.iter_each_modulus(),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, key, ntt_table, modulus)| {
            let (a, b) = glwe.split_at_mut(a_b_mid);

            ntt_table.transform_slice(b);

            let mut res = NttPolynomial(ArrayBase(b));

            a.chunks_exact(poly_length)
                .zip(key.chunks_exact(poly_length))
                .for_each(|(a, s)| {
                    res.add_mul_assign(
                        &NttPolynomial(ArrayBase(a)),
                        &NttPolynomial(ArrayBase(s)),
                        *modulus,
                    );
                });
        });

        Self {
            key: DcrtGlwe::new(ArrayBase(data)),
            moduli_count,
            poly_length,
            dimension,
            a_b_mid,
            glwe_len,
            crt_glwe_len: key_len,
            moduli: moduli.to_vec(),
            modulus_values,
        }
    }

    pub fn encrypt<Table, R, A>(
        &self,
        message: &CrtPolynomial<A>,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        rng: &mut R,
    ) -> DcrtGlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        A: RawData<Elem = T> + Data,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let moduli_count = self.moduli_count;
        let poly_length = self.poly_length;
        let a_b_mid = self.a_b_mid;
        let glwe_len = self.glwe_len;
        let crt_glwe_len = self.crt_glwe_len;

        let mut result = vec![T::ZERO; crt_glwe_len];
        let mut temp = vec![T::ZERO; moduli_count * poly_length];

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut result,
            glwe_len,
            self.modulus_values(),
            gaussian,
            rng,
        );
        primus_distr::sample_crt_binary_values_inplace(&mut temp, poly_length, rng);

        izip!(
            result.chunks_exact_mut(glwe_len),
            self.iter_each_modulus(),
            message.iter_each_modulus(poly_length),
            temp.chunks_exact_mut(poly_length),
            table.iter(),
            self.moduli()
        )
        .for_each(|(glwe, key, msg, v, ntt_table, modulus)| {
            ntt_table.transform_slice(v);
            let v_poly = NttPolynomial(ArrayBase(v));

            Polynomial(ArrayBase(&mut glwe[a_b_mid..]))
                .add_assign(&Polynomial(ArrayBase(msg)), *modulus);

            izip!(
                glwe.chunks_exact_mut(poly_length),
                key.chunks_exact(poly_length)
            )
            .for_each(|(a, k)| {
                ntt_table.transform_slice(a);
                NttPolynomial(ArrayBase(a)).add_mul_assign(
                    &NttPolynomial(ArrayBase(k)),
                    &v_poly,
                    *modulus,
                );
            });
        });

        DcrtGlwe::new(ArrayBase(result))
    }

    fn encrypt_monomial_inner<Table, R, A>(
        &self,
        coeff_residue: &[T],
        degree: usize,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        v: &mut DcrtPolynomial<A>,
        rng: &mut R,
    ) -> DcrtGlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let a_b_mid = self.a_b_mid;
        let glwe_len = self.glwe_len;
        let crt_glwe_len = self.crt_glwe_len;

        let mut result = vec![T::ZERO; crt_glwe_len];

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut result,
            glwe_len,
            self.modulus_values(),
            gaussian,
            rng,
        );
        primus_distr::sample_crt_binary_values_inplace(v.0.as_mut(), poly_length, rng);

        izip!(
            result.chunks_exact_mut(glwe_len),
            self.iter_each_modulus(),
            v.iter_each_modulus_mut(poly_length),
            coeff_residue,
            table.iter(),
            self.moduli()
        )
        .for_each(|(glwe, key, v_r, &coeff, ntt_table, modulus)| {
            ntt_table.transform_slice(v_r);
            let v_poly = NttPolynomial(ArrayBase(v_r));

            glwe[a_b_mid + degree].add_modulo(coeff, *modulus);

            izip!(
                glwe.chunks_exact_mut(poly_length),
                key.chunks_exact(poly_length)
            )
            .for_each(|(a, k)| {
                ntt_table.transform_slice(a);
                NttPolynomial(ArrayBase(a)).add_mul_assign(
                    &NttPolynomial(ArrayBase(k)),
                    &v_poly,
                    *modulus,
                );
            });
        });

        DcrtGlwe::new(ArrayBase(result))
    }

    fn encrypt_neg_secret_monomial_inner<Table, R, A>(
        &self,
        s_index: usize,
        coeff_residue: &[T],
        degree: usize,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        v: &mut DcrtPolynomial<A>,
        rng: &mut R,
    ) -> DcrtGlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let glwe_len = self.glwe_len;
        let crt_glwe_len = self.crt_glwe_len;

        let mut result = vec![T::ZERO; crt_glwe_len];

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut result,
            glwe_len,
            self.modulus_values(),
            gaussian,
            rng,
        );
        primus_distr::sample_crt_binary_values_inplace(v.0.as_mut(), poly_length, rng);

        izip!(
            result.chunks_exact_mut(glwe_len),
            self.iter_each_modulus(),
            v.iter_each_modulus_mut(poly_length),
            coeff_residue,
            table.iter(),
            self.moduli()
        )
        .for_each(|(glwe, key, v_r, &coeff, ntt_table, modulus)| {
            ntt_table.transform_slice(v_r);
            let v_poly = NttPolynomial(ArrayBase(v_r));

            glwe[s_index * poly_length + degree].add_modulo(coeff, *modulus);

            izip!(
                glwe.chunks_exact_mut(poly_length),
                key.chunks_exact(poly_length)
            )
            .for_each(|(a, k)| {
                ntt_table.transform_slice(a);
                NttPolynomial(ArrayBase(a)).add_mul_assign(
                    &NttPolynomial(ArrayBase(k)),
                    &v_poly,
                    *modulus,
                );
            });
        });

        DcrtGlwe::new(ArrayBase(result))
    }

    /// Generate a [`DcrtGgsw`] ciphertext which encrypted `coeff*X^degree`.
    pub fn encrypt_monomial_ggsw<Table, R>(
        &self,
        coeff_residues: &[T],
        degree: usize,
        basis: &BigUintApproxSignedBasis<T>,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        rng: &mut R,
    ) -> DcrtGgsw<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let moduli_count = self.moduli_count;
        let dimension = self.dimension;
        let poly_length = self.poly_length;
        let glwe_len = self.glwe_len;

        let decompose_length = basis.decompose_length();
        let glev_len = decompose_length * glwe_len;
        let ggsw_len = (dimension + 1) * glev_len;

        let v_glev_len = decompose_length * poly_length;
        let v_ggsw_len = (dimension + 1) * v_glev_len;

        let mut dcrt_ggsw = vec![T::ZERO; ggsw_len * moduli_count];

        let mut v_crt_ggsw: Vec<T> =
            primus_distr::sample_crt_binary_values(v_ggsw_len, moduli_count, rng);

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut dcrt_ggsw,
            ggsw_len,
            self.modulus_values(),
            gaussian,
            rng,
        );

        izip!(
            dcrt_ggsw.chunks_exact_mut(ggsw_len),
            self.iter_each_modulus(),
            v_crt_ggsw.chunks_exact_mut(v_ggsw_len),
            coeff_residues,
            basis.scalars_residue().chunks_exact(decompose_length),
            table.iter(),
            self.moduli()
        )
        .for_each(|(ggsw, key, v_ggsw, &coeff, scalars, ntt_table, modulus)| {
            ggsw.chunks_exact_mut(glev_len)
                .zip(v_ggsw.chunks_exact_mut(v_glev_len))
                .enumerate()
                .for_each(|(i, (glev, v_glev))| {
                    izip!(
                        glev.chunks_exact_mut(glwe_len),
                        v_glev.chunks_exact_mut(poly_length),
                        scalars
                    )
                    .for_each(|(glwe, v_glwe, &scalar)| {
                        let index = i * poly_length + degree;
                        glwe[index] = coeff.mul_add_modulo(scalar, glwe[index], *modulus);

                        ntt_table.transform_slice(v_glwe);
                        let v_poly = NttPolynomial(ArrayBase(&*v_glwe));

                        glwe.chunks_exact_mut(poly_length)
                            .zip(key.chunks_exact(poly_length))
                            .for_each(|(a, s)| {
                                ntt_table.transform_slice(a);
                                NttPolynomial(ArrayBase(a)).add_mul_assign(
                                    &v_poly,
                                    &NttPolynomial(ArrayBase(s)),
                                    *modulus,
                                );
                            });
                    });
                });
        });

        DcrtGgsw {
            data: ArrayBase(dcrt_ggsw),
        }
    }
}
