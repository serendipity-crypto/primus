use std::sync::Arc;

use primus_factor::ShoupFactor;
use primus_integer::{AsInto, UnsignedInteger};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, RawData};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{
    CrtGlevParameters, CrtGlweAutoContext, CrtGlweAutoKey, CrtGlweCiphertext, CrtGlweSecretKey,
    DcrtGlweCiphertext, DcrtGlweSecretKey,
};

pub struct CrtGlweTraceContext<T: UnsignedInteger> {
    crt_glwe: CrtGlweCiphertext<Vec<T>>,
    auto_context: CrtGlweAutoContext<T>,
}

impl<T: UnsignedInteger> CrtGlweTraceContext<T> {
    pub fn new(
        dimension: usize,
        poly_length: usize,
        crt_poly_len: usize,
        big_uint_poly_len: usize,
    ) -> Self {
        let crt_glwe = CrtGlweCiphertext::zero((dimension + 1) * crt_poly_len);
        let auto_context = CrtGlweAutoContext::new(poly_length, crt_poly_len, big_uint_poly_len);
        Self {
            crt_glwe,
            auto_context,
        }
    }

    pub fn as_mut(
        &mut self,
    ) -> (
        &mut primus_lattice::glwe::CrtGlwe<Vec<T>, T>,
        &mut CrtGlweAutoContext<T>,
    ) {
        (&mut self.crt_glwe, &mut self.auto_context)
    }
}

#[derive(Clone)]
pub struct CrtGlweTraceKey<T: UnsignedInteger, Table>
where
    Table: DcrtTable<ValueT = T> + Dcrt,
{
    auto_keys: Vec<CrtGlweAutoKey<T, Table>>,
    table: Arc<Table>,
}

impl<T: UnsignedInteger, Table> CrtGlweTraceKey<T, Table>
where
    Table: DcrtTable<ValueT = T> + Dcrt,
{
    pub fn new<M, R>(
        params: &CrtGlevParameters<T, M>,
        sk: &CrtGlweSecretKey<T>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        table: Arc<Table>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let log_n = params.poly_length().trailing_zeros();
        let auto_keys: Vec<CrtGlweAutoKey<T, Table>> = (1..=log_n)
            .rev()
            .map(|x| (1usize << x) + 1)
            .map(|degree| CrtGlweAutoKey::new(params, degree, sk, dcrt_sk, Arc::clone(&table), rng))
            .collect();
        Self { auto_keys, table }
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    pub fn trace_inplace<M, A, B>(
        &self,
        ciphertext: &CrtGlweCiphertext<A>,
        result: &mut CrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweTraceContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let moduli = params.cipher_moduli();

        let (crt_glwe, auto_context) = context.as_mut();

        result.as_mut().copy_from_slice(ciphertext.as_ref());

        for auto_key in self.auto_keys.iter() {
            auto_key.automorphism_inplace(result, crt_glwe, params, rns_base, auto_context);
            result.add_element_wise_assign(crt_glwe, poly_length, rns_poly_len, moduli);
        }
    }

    pub fn trace_to_dcrt_glwe_inplace<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut DcrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweTraceContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let moduli = params.cipher_moduli();

        let (crt_glwe, auto_context) = context.as_mut();

        let dcrt_glwe = &mut DcrtGlweCiphertext::new(ArrayBase(crt_glwe.as_mut()));

        result.as_mut().copy_from_slice(ciphertext.as_ref());

        for auto_key in self.auto_keys.iter() {
            auto_key.automorphism_to_dcrt_glwe_inplace(
                result,
                dcrt_glwe,
                params,
                rns_base,
                auto_context,
            );
            result.add_element_wise_assign(dcrt_glwe, poly_length, rns_poly_len, moduli);
        }
    }

    /// Coefficient Expansion Algorithm.
    ///
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_coefficients_inplace<M, A, B>(
        &self,
        ciphertext: &CrtGlweCiphertext<A>,
        result: &mut [CrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweTraceContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let twice_poly_length = poly_length * 2;
        let rns_poly_len = params.rns_poly_len();
        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();

        let big_uint_value_len = params.big_uint_value_len();
        let mut n = vec![T::ZERO; big_uint_value_len];
        n[0] = poly_length.as_into();
        let n_residue = rns_base.decompose(&n);

        let inv_n_residue: Vec<ShoupFactor<T>> = n_residue
            .iter()
            .zip(rns_base.moduli())
            .map(|(&n, m)| ShoupFactor::new(m.reduce_inv(n), m.value_unchecked()))
            .collect();

        ciphertext.mul_factor_inplace(
            &inv_n_residue,
            &mut result[0],
            poly_length,
            rns_poly_len,
            moduli_value,
        );

        let (crt_glwe, auto_context) = context.as_mut();

        for (i, auto_key) in self.auto_keys.iter().enumerate() {
            let two_pow_i = 1 << i;

            let (x, y) = unsafe { result[..two_pow_i * 2].split_at_mut_unchecked(two_pow_i) };

            x.iter_mut().zip(y.iter_mut()).for_each(|(a_0, b_0)| {
                auto_key.automorphism_inplace(a_0, crt_glwe, params, rns_base, auto_context);

                a_0.sub_element_wise_inplace(crt_glwe, b_0, poly_length, rns_poly_len, moduli);
                b_0.mul_monic_monomial_assign(
                    twice_poly_length - two_pow_i,
                    poly_length,
                    rns_poly_len,
                    moduli,
                );
                a_0.add_element_wise_assign(crt_glwe, poly_length, rns_poly_len, moduli);
            });
        }
    }

    /// Coefficient Expansion Algorithm.
    ///
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_partial_coefficients_inplace<M, A, B>(
        &self,
        ciphertext: &CrtGlweCiphertext<A>,
        result: &mut [CrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweTraceContext<T>,
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
        let big_uint_value_len = params.big_uint_value_len();

        let twice_poly_length = poly_length * 2;

        let log_d = count.trailing_zeros() as usize;
        let mut n = vec![T::ZERO; big_uint_value_len];
        n[0] = count.as_into();
        let n_residue = rns_base.decompose(&n);

        let inv_n_residue: Vec<ShoupFactor<T>> = n_residue
            .iter()
            .zip(rns_base.moduli())
            .map(|(&n, m)| ShoupFactor::new(m.reduce_inv(n), m.value_unchecked()))
            .collect();

        ciphertext.mul_factor_inplace(
            &inv_n_residue,
            &mut result[0],
            poly_length,
            rns_poly_len,
            moduli_value,
        );

        let (crt_glwe, auto_context) = context.as_mut();

        for (i, auto_key) in self.auto_keys.iter().enumerate().take(log_d) {
            let two_pow_i = 1 << i;

            let (x, y) = unsafe { result[..two_pow_i * 2].split_at_mut_unchecked(two_pow_i) };

            x.iter_mut().zip(y.iter_mut()).for_each(|(a_0, b_0)| {
                auto_key.automorphism_inplace(a_0, crt_glwe, params, rns_base, auto_context);

                a_0.sub_element_wise_inplace(crt_glwe, b_0, poly_length, rns_poly_len, moduli);
                b_0.mul_monic_monomial_assign(
                    twice_poly_length - two_pow_i,
                    poly_length,
                    rns_poly_len,
                    moduli,
                );
                a_0.add_element_wise_assign(crt_glwe, poly_length, rns_poly_len, moduli);
            });
        }
    }
}
