use std::sync::Arc;

use primus_factor::ShoupFactor;
use primus_integer::{AsInto, BigUint, Data, DataMut, RawData, UnsignedInteger};
use primus_ntt::DcrtTable;
use primus_poly::DcrtPolynomial;
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{
    CrtGlevParameters, DcrtGlweAutoKey, DcrtGlweCiphertext, DcrtGlweSecretKey, DcrtGlweTraceContext,
};

pub type DcrtGlweExpandCoeffContext<T> = DcrtGlweTraceContext<T>;

pub struct DcrtGlweExpandCoeffContextPool<T: UnsignedInteger> {
    contexts: Vec<DcrtGlweExpandCoeffContext<T>>,
}

impl<T: UnsignedInteger> DcrtGlweExpandCoeffContextPool<T> {
    pub fn new(
        size: usize,
        dimension: usize,
        poly_length: usize,
        crt_poly_len: usize,
        big_uint_poly_len: usize,
    ) -> Self {
        assert!(size > 0);
        let contexts = (0..size)
            .map(|_| {
                DcrtGlweExpandCoeffContext::new(
                    dimension,
                    poly_length,
                    crt_poly_len,
                    big_uint_poly_len,
                )
            })
            .collect();
        Self { contexts }
    }

    pub fn get_mut(&mut self, idx: usize) -> &mut DcrtGlweExpandCoeffContext<T> {
        &mut self.contexts[idx]
    }

    pub fn as_mut_slice(&mut self) -> &mut [DcrtGlweExpandCoeffContext<T>] {
        self.contexts.as_mut_slice()
    }

    pub fn len(&self) -> usize {
        self.contexts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contexts.is_empty()
    }
}

#[derive(Clone)]
pub struct DcrtGlweExpandCoeffKey<T: UnsignedInteger, Table>
where
    Table: DcrtTable<ValueT = T>,
{
    auto_keys: Vec<DcrtGlweAutoKey<T, Table>>,
    monomial_ntt_by_level: Vec<Vec<T>>,
    inv_n_residue_by_log_count: Vec<Vec<ShoupFactor<T>>>,
    table: Arc<Table>,
}

impl<T: UnsignedInteger, Table> DcrtGlweExpandCoeffKey<T, Table>
where
    Table: DcrtTable<ValueT = T>,
{
    pub fn new<M, R>(
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        table: Arc<Table>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let log_n = params.poly_length().trailing_zeros();
        let auto_keys: Vec<DcrtGlweAutoKey<T, Table>> = (1..=log_n)
            .rev()
            .map(|x| (1usize << x) + 1)
            .map(|degree| DcrtGlweAutoKey::new(params, degree, dcrt_sk, Arc::clone(&table), rng))
            .collect();

        let monomial_ntt_by_level = Self::precompute_monomials(params, table.as_ref());
        let inv_n_residue_by_log_count = Self::precompute_inv_n_residue(params, rns_base);

        Self {
            auto_keys,
            monomial_ntt_by_level,
            inv_n_residue_by_log_count,
            table,
        }
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    fn precompute_monomials<M>(params: &CrtGlevParameters<T, M>, table: &Table) -> Vec<Vec<T>>
    where
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let twice_poly_length = poly_length * 2;
        let log_n = poly_length.trailing_zeros() as usize;

        (0..log_n)
            .map(|i| {
                let degree = twice_poly_length - (1 << i);
                let mut monomial_ntt = vec![T::ZERO; rns_poly_len];
                table.transform_coeff_one_monomial(degree, &mut monomial_ntt);
                monomial_ntt
            })
            .collect()
    }

    fn precompute_inv_n_residue<M>(
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
    ) -> Vec<Vec<ShoupFactor<T>>>
    where
        M: FieldContext<T>,
    {
        let big_uint_value_len = params.big_uint_value_len();
        let log_n = params.poly_length().trailing_zeros() as usize;

        (0..=log_n)
            .map(|log_count| {
                let count = 1usize << log_count;
                let mut n = vec![T::ZERO; big_uint_value_len];
                n[0] = count.as_into();
                let n_residue = rns_base.decompose(BigUint(&n));

                n_residue
                    .iter()
                    .zip(rns_base.moduli())
                    .map(|(&n, m)| ShoupFactor::new(m.reduce_inv(n), m.value_unchecked()))
                    .collect()
            })
            .collect()
    }

    /// Coefficient Expansion Algorithm.
    ///
    /// Expands all `poly_length` coefficients.
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_coefficients_inplace<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut [DcrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut DcrtGlweExpandCoeffContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        debug_assert_eq!(result.len(), params.poly_length());
        self.expand_partial_coefficients_inplace(ciphertext, result, params, rns_base, context);
    }

    /// Coefficient Expansion Algorithm.
    ///
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_partial_coefficients_inplace<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut [DcrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut DcrtGlweExpandCoeffContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let count = result.len();
        assert!(count.is_power_of_two() && count <= poly_length);

        let rns_poly_len = params.rns_poly_len();
        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();

        let log_d = count.trailing_zeros() as usize;
        debug_assert!(log_d <= self.monomial_ntt_by_level.len());
        debug_assert!(log_d < self.inv_n_residue_by_log_count.len());

        ciphertext.mul_factor_inplace(
            &self.inv_n_residue_by_log_count[log_d],
            &mut result[0],
            poly_length,
            rns_poly_len,
            moduli_value,
        );

        let (dcrt_glwe, auto_context) = context.as_mut();

        for (i, auto_key) in self.auto_keys.iter().enumerate().take(log_d) {
            let two_pow_i = 1 << i;
            let monomial_ntt_poly = DcrtPolynomial(self.monomial_ntt_by_level[i].as_slice());

            let (x, y) = unsafe { result[..two_pow_i * 2].split_at_mut_unchecked(two_pow_i) };

            x.iter_mut().zip(y.iter_mut()).for_each(|(a_0, b_0)| {
                auto_key.automorphism_inplace(a_0, dcrt_glwe, params, rns_base, auto_context);

                a_0.butterfly_mul_dcrt_polynomial_inplace(
                    dcrt_glwe,
                    &monomial_ntt_poly,
                    b_0,
                    poly_length,
                    moduli,
                );
            });
        }
    }
}
