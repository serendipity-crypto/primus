use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::{ggsw::DcrtGgsw, glwe::DcrtGlwe};
use primus_modulo::AddModulo;
use primus_modulo::MulAddModulo;
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, NttPolynomial, Polynomial, RawData, crt::CrtPolynomial,
    dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::CrtGlweParameters;
use crate::DcrtGlweSecretKey;

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

        let key_len = crt_poly_length * dimension;

        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let mut data = vec![T::ZERO; key_len];

        let (a, b) = unsafe { data.split_at_mut_unchecked(key_len) };

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
            key: DcrtGlwe::new(ArrayBase(data)),
            moduli_count: secret_key.moduli_count(),
            poly_length,
            dimension,
            crt_poly_length,
            crt_glwe_len: crt_poly_length * (dimension + 1),
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
        let moduli_count = self.moduli_count;
        let poly_length = self.poly_length;
        todo!()
        // let a_b_mid = self.a_b_mid;
        // let glwe_len = self.glwe_len;
        // let crt_glwe_len = self.crt_glwe_len;

        // let mut result = vec![T::ZERO; crt_glwe_len];
        // let mut temp = vec![T::ZERO; moduli_count * poly_length];

        // primus_distr::sample_crt_gaussian_values_inplace(
        //     &mut result,
        //     glwe_len,
        //     self.modulus_values(),
        //     gaussian,
        //     rng,
        // );
        // primus_distr::sample_crt_binary_values_inplace(&mut temp, poly_length, rng);

        // izip!(
        //     result.chunks_exact_mut(glwe_len),
        //     self.iter_each_modulus(),
        //     message.iter_each_modulus(poly_length),
        //     temp.chunks_exact_mut(poly_length),
        //     table.iter(),
        //     self.moduli()
        // )
        // .for_each(|(glwe, key, msg, v, ntt_table, modulus)| {
        //     ntt_table.transform_slice(v);
        //     let v_poly = NttPolynomial(ArrayBase(v));

        //     Polynomial(ArrayBase(&mut glwe[a_b_mid..]))
        //         .add_assign(&Polynomial(ArrayBase(msg)), *modulus);

        //     izip!(
        //         glwe.chunks_exact_mut(poly_length),
        //         key.chunks_exact(poly_length)
        //     )
        //     .for_each(|(a, k)| {
        //         ntt_table.transform_slice(a);
        //         NttPolynomial(ArrayBase(a)).add_mul_assign(
        //             &NttPolynomial(ArrayBase(k)),
        //             &v_poly,
        //             *modulus,
        //         );
        //     });
        // });

        // DcrtGlwe::new(ArrayBase(result))
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
        todo!()
        // let a_b_mid = self.a_b_mid;
        // let glwe_len = self.glwe_len;
        // let crt_glwe_len = self.crt_glwe_len;

        // let mut result = vec![T::ZERO; crt_glwe_len];

        // primus_distr::sample_crt_gaussian_values_inplace(
        //     &mut result,
        //     glwe_len,
        //     self.modulus_values(),
        //     gaussian,
        //     rng,
        // );
        // primus_distr::sample_crt_binary_values_inplace(v.0.as_mut(), poly_length, rng);

        // izip!(
        //     result.chunks_exact_mut(glwe_len),
        //     self.iter_each_modulus(),
        //     v.iter_each_modulus_mut(poly_length),
        //     coeff_residue,
        //     table.iter(),
        //     self.moduli()
        // )
        // .for_each(|(glwe, key, v_r, &coeff, ntt_table, modulus)| {
        //     ntt_table.transform_slice(v_r);
        //     let v_poly = NttPolynomial(ArrayBase(v_r));

        //     glwe[a_b_mid + degree].add_modulo(coeff, *modulus);

        //     izip!(
        //         glwe.chunks_exact_mut(poly_length),
        //         key.chunks_exact(poly_length)
        //     )
        //     .for_each(|(a, k)| {
        //         ntt_table.transform_slice(a);
        //         NttPolynomial(ArrayBase(a)).add_mul_assign(
        //             &NttPolynomial(ArrayBase(k)),
        //             &v_poly,
        //             *modulus,
        //         );
        //     });
        // });

        // DcrtGlwe::new(ArrayBase(result))
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
        todo!()
        // let glwe_len = self.glwe_len;
        // let crt_glwe_len = self.crt_glwe_len;

        // let mut result = vec![T::ZERO; crt_glwe_len];

        // primus_distr::sample_crt_gaussian_values_inplace(
        //     &mut result,
        //     glwe_len,
        //     self.modulus_values(),
        //     gaussian,
        //     rng,
        // );
        // primus_distr::sample_crt_binary_values_inplace(v.0.as_mut(), poly_length, rng);

        // izip!(
        //     result.chunks_exact_mut(glwe_len),
        //     self.iter_each_modulus(),
        //     v.iter_each_modulus_mut(poly_length),
        //     coeff_residue,
        //     table.iter(),
        //     self.moduli()
        // )
        // .for_each(|(glwe, key, v_r, &coeff, ntt_table, modulus)| {
        //     ntt_table.transform_slice(v_r);
        //     let v_poly = NttPolynomial(ArrayBase(v_r));

        //     glwe[s_index * poly_length + degree].add_modulo(coeff, *modulus);

        //     izip!(
        //         glwe.chunks_exact_mut(poly_length),
        //         key.chunks_exact(poly_length)
        //     )
        //     .for_each(|(a, k)| {
        //         ntt_table.transform_slice(a);
        //         NttPolynomial(ArrayBase(a)).add_mul_assign(
        //             &NttPolynomial(ArrayBase(k)),
        //             &v_poly,
        //             *modulus,
        //         );
        //     });
        // });

        // DcrtGlwe::new(ArrayBase(result))
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
        todo!()
        // let dimension = self.dimension;
        // let poly_length = self.poly_length;
        // let glwe_len = self.glwe_len;

        // let decompose_length = basis.decompose_length();
        // let glev_len = decompose_length * glwe_len;
        // let ggsw_len = (dimension + 1) * glev_len;

        // let v_glev_len = decompose_length * poly_length;
        // let v_ggsw_len = (dimension + 1) * v_glev_len;

        // let mut dcrt_ggsw = vec![T::ZERO; ggsw_len * moduli_count];

        // let mut v_crt_ggsw: Vec<T> =
        //     primus_distr::sample_crt_binary_values(v_ggsw_len, moduli_count, rng);

        // primus_distr::sample_crt_gaussian_values_inplace(
        //     &mut dcrt_ggsw,
        //     ggsw_len,
        //     self.modulus_values(),
        //     gaussian,
        //     rng,
        // );

        // izip!(
        //     dcrt_ggsw.chunks_exact_mut(ggsw_len),
        //     self.iter_each_modulus(),
        //     v_crt_ggsw.chunks_exact_mut(v_ggsw_len),
        //     coeff_residues,
        //     basis.scalars_residue().chunks_exact(decompose_length),
        //     table.iter(),
        //     self.moduli()
        // )
        // .for_each(|(ggsw, key, v_ggsw, &coeff, scalars, ntt_table, modulus)| {
        //     ggsw.chunks_exact_mut(glev_len)
        //         .zip(v_ggsw.chunks_exact_mut(v_glev_len))
        //         .enumerate()
        //         .for_each(|(i, (glev, v_glev))| {
        //             izip!(
        //                 glev.chunks_exact_mut(glwe_len),
        //                 v_glev.chunks_exact_mut(poly_length),
        //                 scalars
        //             )
        //             .for_each(|(glwe, v_glwe, &scalar)| {
        //                 let index = i * poly_length + degree;
        //                 glwe[index] = coeff.mul_add_modulo(scalar, glwe[index], *modulus);

        //                 ntt_table.transform_slice(v_glwe);
        //                 let v_poly = NttPolynomial(ArrayBase(&*v_glwe));

        //                 glwe.chunks_exact_mut(poly_length)
        //                     .zip(key.chunks_exact(poly_length))
        //                     .for_each(|(a, s)| {
        //                         ntt_table.transform_slice(a);
        //                         NttPolynomial(ArrayBase(a)).add_mul_assign(
        //                             &v_poly,
        //                             &NttPolynomial(ArrayBase(s)),
        //                             *modulus,
        //                         );
        //                     });
        //             });
        //         });
        // });

        // DcrtGgsw {
        //     data: ArrayBase(dcrt_ggsw),
        // }
    }
}
