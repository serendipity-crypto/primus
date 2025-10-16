use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_reduce::FieldContext;

use crate::{CrtGlweSecretKey, DcrtGlweSecretKey};

pub struct DcrtGlweKeySwitchingKey<T: UnsignedInteger> {
    key: Vec<T>,
}

impl<T: UnsignedInteger> DcrtGlweKeySwitchingKey<T> {
    pub fn new_auto_key<R, M, Table>(
        sk: &CrtGlweSecretKey<T>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        degree: usize,
        gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
        moduli: &[M],
        table: &Table,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = sk.poly_length();
        let dimension = sk.dimension();
        let moduli_count = sk.moduli_count();

        let a_b_mid = dimension * poly_length;
        let glwe_len = (dimension + 1) * poly_length;
        let single_modulus_len = dimension * glwe_len;

        let e_single_modulus_len = a_b_mid;

        let mut result = vec![T::ZERO; moduli_count * single_modulus_len];
        let mut e_all = vec![T::ZERO; moduli_count * e_single_modulus_len];
        let modulus_values: Vec<T> = moduli.iter().map(|m| m.value_unchecked()).collect();

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut e_all,
            e_single_modulus_len,
            &modulus_values,
            gaussian,
            rng,
        );

        izip!(
            result.chunks_exact_mut(single_modulus_len),
            sk.iter_each_modulus(),
            e_all.chunks_exact_mut(e_single_modulus_len),
            table.iter(),
            moduli
        )
        .for_each(|(auto_key, key, es, ntt_table, modulus)| {
            izip!(
                auto_key.chunks_exact_mut(glwe_len),
                es.chunks_exact_mut(poly_length)
            )
            .for_each(|(glwe, e)| {
                let (a, b) = unsafe { glwe.split_at_mut_unchecked(e_single_modulus_len) };
                ntt_table.transform_slice(e);
            });
        });
        todo!()
    }
}
