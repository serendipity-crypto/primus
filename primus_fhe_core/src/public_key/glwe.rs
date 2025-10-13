use std::slice::IterMut;

use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::DcrtGlwe;
use primus_modulo::AddModulo;
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

        let mid = poly_length * dimension;
        let glwe_len = mid + poly_length;
        let key_len = moduli_count * poly_length * (dimension + 1);

        let modulus_values: Vec<T> = moduli.iter().map(|m| m.value_unchecked()).collect();
        let uniform_distrs: Vec<Uniform<T>> =
            moduli.iter().map(|m| m.uniform_distribution()).collect();

        let mut data = vec![T::ZERO; key_len];

        let (mut a, mut b): (Vec<IterMut<'_, T>>, Vec<IterMut<'_, T>>) = data
            .chunks_exact_mut(glwe_len)
            .map(|s| {
                let (a, b) = s.split_at_mut(mid);
                (a.iter_mut(), b.iter_mut())
            })
            .collect();

        primus_distr::sample_crt_uniform_values_iter_mut(&mut a, &uniform_distrs, rng);
        primus_distr::sample_crt_gaussian_values_iter_mut(&mut b, &modulus_values, gaussian, rng);

        izip!(data.chunks_exact_mut(glwe_len), table.iter(), moduli).for_each(
            |(glwe, ntt_table, modulus)| {
                let (a, b) = glwe.split_at_mut(mid);

                ntt_table.transform_slice(b);

                let mut res = NttPolynomial(ArrayBase(b));

                a.chunks_exact(poly_length)
                    .zip(secret_key.key.chunks_exact(poly_length))
                    .for_each(|(a, s)| {
                        res.add_mul_assign(
                            &NttPolynomial(ArrayBase(a)),
                            &NttPolynomial(ArrayBase(s)),
                            *modulus,
                        );
                    });
            },
        );

        Self {
            key: DcrtGlwe::new(ArrayBase(data)),
            moduli_count,
            poly_length,
            dimension,
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
        let dimension = self.dimension;

        let mid = poly_length * dimension;
        let glwe_len = mid + poly_length;
        let key_len = moduli_count * glwe_len;

        let mut result = vec![T::ZERO; key_len];
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
            self.key.data.chunks_exact(glwe_len),
            message.0.chunks_exact(poly_length),
            temp.chunks_exact_mut(poly_length),
            table.iter(),
            self.moduli()
        )
        .for_each(|(glwe, key, msg, v, ntt_table, modulus)| {
            ntt_table.transform_slice(v);
            let v_poly = NttPolynomial(ArrayBase(v));

            Polynomial(ArrayBase(&mut glwe[mid..]))
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
        let moduli_count = self.moduli_count;
        let poly_length = self.poly_length;
        let dimension = self.dimension;

        let mid = poly_length * dimension;
        let glwe_len = mid + poly_length;
        let key_len = moduli_count * glwe_len;

        let mut result = vec![T::ZERO; key_len];
        // let mut temp = vec![T::ZERO; moduli_count * poly_length];

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
            self.key.data.chunks_exact(glwe_len),
            v.iter_mut(poly_length),
            coeff_residue,
            table.iter(),
            self.moduli()
        )
        .for_each(|(glwe, key, v_r, &residue, ntt_table, modulus)| {
            ntt_table.transform_slice(v_r);
            let v_poly = NttPolynomial(ArrayBase(v_r));

            glwe[mid + degree].add_modulo(residue, *modulus);

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
}
