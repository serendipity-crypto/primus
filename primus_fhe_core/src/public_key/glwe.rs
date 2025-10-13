use std::slice::IterMut;

use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::DcrtGlwe;
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, NttPolynomial, RawData, crt::CrtPolynomial};
use primus_reduce::FieldContext;
use rand::distr::Uniform;

use crate::DcrtGlweSecretKey;

pub struct DcrtGlwePublicKey<T: UnsignedInteger> {
    key: DcrtGlwe<Vec<T>>,
    moduli_count: usize,
    poly_length: usize,
    dimension: usize,
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

    pub fn new<Table, M, R>(
        secret_key: &DcrtGlweSecretKey<T>,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        moduli: &[M],
        rng: &mut R,
    ) -> DcrtGlwePublicKey<T>
    where
        R: rand::Rng + rand::CryptoRng,
        Table: DcrtTable<ValueT = T> + Dcrt,
        M: FieldContext<T>,
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
        }
    }

    pub fn encrypt<Table, M, R, A>(
        &self,
        message: &CrtPolynomial<A>,
        gaussian: &SignedDiscreteGaussian<T::SignedInteger>,
        table: &Table,
        moduli: &[M],
        rng: &mut R,
    ) -> DcrtGlwe<Vec<T>>
    where
        R: rand::Rng + rand::CryptoRng,
        A: RawData<Elem = T> + Data,
        Table: DcrtTable<ValueT = T> + Dcrt,
        M: FieldContext<T>,
    {
        let moduli_count = self.moduli_count;
        let poly_length = self.poly_length;
        let dimension = self.dimension;

        let mid = poly_length * dimension;
        let glwe_len = mid + poly_length;
        let key_len = moduli_count * glwe_len;

        let mut result = vec![T::ZERO; key_len];
        let mut temp = vec![T::ZERO; moduli_count * poly_length];
        let modulus_values: Vec<T> = moduli.iter().map(|m| m.value_unchecked()).collect();

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut result,
            glwe_len,
            &modulus_values,
            gaussian,
            rng,
        );
        primus_distr::sample_crt_binary_values_inplace(&mut temp, poly_length, rng);

        izip!(
            result.chunks_exact_mut(glwe_len),
            self.key.data.chunks_exact(glwe_len),
            temp.chunks_exact_mut(poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, key, v, ntt_table, modulus)| {
            ntt_table.transform_slice(v);
            let v_poly = NttPolynomial(ArrayBase(v));
            izip!(
                glwe.chunks_exact_mut(poly_length),
                key.chunks_exact(poly_length)
            )
            .for_each(|(a, k)| {
                NttPolynomial(ArrayBase(a)).add_mul_assign(
                    &NttPolynomial(ArrayBase(k)),
                    &v_poly,
                    *modulus,
                );
            });
        });

        temp.copy_from_slice(message.0.as_ref());

        izip!(
            result.chunks_exact_mut(glwe_len),
            temp.chunks_exact_mut(poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, v, ntt_table, modulus)| {
            ntt_table.transform_slice(v);
            let v_poly = NttPolynomial(ArrayBase(v));
            glwe.chunks_exact_mut(poly_length).for_each(|a| {
                NttPolynomial(ArrayBase(a)).add_assign(&v_poly, *modulus);
            });
        });

        DcrtGlwe::new(ArrayBase(result))
    }
}
