use std::sync::{Arc, Mutex};

use primus_factor::ShoupFactor;
use primus_integer::{AsInto, BigUint, Data, DataMut, RawData, UnsignedInteger};
use primus_ntt::DcrtTable;
use primus_poly::DcrtPolynomial;
use primus_reduce::FieldContext;
use primus_rns::RNSBase;
use rayon::prelude::*;

use crate::{
    CrtGlevParameters, DcrtGlweAutoKey, DcrtGlweCiphertext, DcrtGlweSecretKey, DcrtGlweTraceContext,
};

/// Wrapper that asserts `Sync` for shared immutable references in parallel contexts.
///
/// # Safety
///
/// The wrapped reference must only be used for immutable access. The non-`Sync` fields
/// (e.g., `rand::Uniform` samplers within `CrtGlevParameters`) must not be accessed
/// concurrently through this wrapper.
struct SyncRef<'a, T: ?Sized>(&'a T);

// Safety: The wrapped reference is only shared immutably across threads.
// `CrtGlevParameters` contains `SignedDiscreteGaussian` whose `rand::Uniform` sampler
// doesn't implement `Sync` due to overly conservative generic bounds in the `rand` crate,
// but the sampler fields are never accessed during the parallel coefficient expansion.
unsafe impl<T: ?Sized> Sync for SyncRef<'_, T> {}
unsafe impl<T: ?Sized> Send for SyncRef<'_, T> {}

impl<T: ?Sized> std::ops::Deref for SyncRef<'_, T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.0
    }
}

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

/// Thread-safe context pool for parallel coefficient expansion.
///
/// Contexts are lazily allocated on demand and reused via `acquire`/`release`.
/// The pool grows up to the number of concurrent worker threads.
pub struct DcrtGlweExpandCoeffSyncPool<T: UnsignedInteger> {
    contexts: Mutex<Vec<DcrtGlweExpandCoeffContext<T>>>,
    dimension: usize,
    poly_length: usize,
    crt_poly_len: usize,
    big_uint_poly_len: usize,
}

impl<T: UnsignedInteger> DcrtGlweExpandCoeffSyncPool<T> {
    /// Creates an empty pool. Contexts are allocated lazily on first [`Self::acquire`].
    pub fn new(
        dimension: usize,
        poly_length: usize,
        crt_poly_len: usize,
        big_uint_poly_len: usize,
    ) -> Self {
        Self {
            contexts: Mutex::new(Vec::new()),
            dimension,
            poly_length,
            crt_poly_len,
            big_uint_poly_len,
        }
    }

    /// Pop a context from the pool, or create a new one if empty.
    pub fn acquire(&self) -> DcrtGlweExpandCoeffContext<T> {
        self.contexts.lock().unwrap().pop().unwrap_or_else(|| {
            DcrtGlweExpandCoeffContext::new(
                self.dimension,
                self.poly_length,
                self.crt_poly_len,
                self.big_uint_poly_len,
            )
        })
    }

    /// Return a context to the pool for reuse.
    pub fn release(&self, ctx: DcrtGlweExpandCoeffContext<T>) {
        self.contexts.lock().unwrap().push(ctx);
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

    /// Parallel Coefficient Expansion Algorithm.
    ///
    /// Expands all `poly_length` coefficients using rayon parallelism.
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_coefficients_inplace_parallel<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut [DcrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context_pool: &DcrtGlweExpandCoeffSyncPool<T>,
    ) where
        M: FieldContext<T> + Sync,
        A: RawData<Elem = T> + Data + Sync,
        B: RawData<Elem = T> + DataMut + Send,
        Table: Send + Sync,
    {
        debug_assert_eq!(result.len(), params.poly_length());
        self.expand_partial_coefficients_inplace_parallel(
            ciphertext,
            result,
            params,
            rns_base,
            context_pool,
        );
    }

    /// Parallel Coefficient Expansion Algorithm.
    ///
    /// (Alg. 1)<https://eprint.iacr.org/2024/266.pdf>
    pub fn expand_partial_coefficients_inplace_parallel<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut [DcrtGlweCiphertext<B>],
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context_pool: &DcrtGlweExpandCoeffSyncPool<T>,
    ) where
        M: FieldContext<T> + Sync,
        A: RawData<Elem = T> + Data + Sync,
        B: RawData<Elem = T> + DataMut + Send,
        Table: Send + Sync,
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

        // Wrap params and rns_base in SyncRef to allow sharing across rayon threads.
        // CrtGlevParameters contains SignedDiscreteGaussian whose rand::Uniform sampler
        // doesn't implement Sync due to generic bounds in the rand crate, but those
        // fields are never accessed during the parallel butterfly computation.
        let params = SyncRef(params);
        let rns_base = SyncRef(rns_base);

        for (i, auto_key) in self.auto_keys.iter().enumerate().take(log_d) {
            let two_pow_i = 1 << i;
            let monomial_ntt_poly = DcrtPolynomial(self.monomial_ntt_by_level[i].as_slice());

            let (x, y) = unsafe { result[..two_pow_i * 2].split_at_mut_unchecked(two_pow_i) };

            x.par_iter_mut()
                .zip(y.par_iter_mut())
                .for_each(|(a_0, b_0)| {
                    let mut ctx = context_pool.acquire();
                    let (dcrt_glwe, auto_context) = ctx.as_mut();

                    auto_key.automorphism_inplace(a_0, dcrt_glwe, &params, &rns_base, auto_context);

                    a_0.butterfly_mul_dcrt_polynomial_inplace(
                        dcrt_glwe,
                        &monomial_ntt_poly,
                        b_0,
                        poly_length,
                        moduli,
                    );

                    context_pool.release(ctx);
                });
        }
    }
}
