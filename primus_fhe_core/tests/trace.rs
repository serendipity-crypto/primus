use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweParameters, CrtGlweSecretKey, CrtGlweTraceContext, CrtGlweTraceKey,
    DcrtGlweCiphertext, DcrtGlweDecryptContext, DcrtGlweSecretKey, DcrtGlweTraceContext,
    DcrtGlweTraceKey, RingSecretKeyType,
};
use primus_lattice::glwe::CrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{CrtPolynomial, Polynomial};
use primus_reduce::ops::*;

#[test]
fn test_crt_glwe_trace() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 1 << 15;
    // let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    // let gamma: ValueT = 2305843009213554689;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; _] = [1125899906826241, 1125899906629633];
    let moduli = moduli_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(log_n, &moduli).unwrap();

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        dimension,
        poly_length,
        mod_t,
        mod_gamma,
        &moduli,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let moduli_count = glwe_params.cipher_moduli_count();
    let rns_poly_len = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let trace_key = CrtGlweTraceKey::new(&glev_params, &sk, &dcrt_sk, Arc::new(table), &mut rng);
    let table = trace_key.table();

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c2: CrtGlwe<Vec<ValueT>> = CrtGlwe::zero(rns_glwe_len);
    let mut trace_context =
        CrtGlweTraceContext::new(dimension, poly_length, rns_poly_len, big_uint_poly_len);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    let mut c1 = c1.into_coeff_form(table);

    trace_key.trace_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut trace_context,
    );

    let c2 = c2.into_ntt_form(table);

    let trace_msg = dcrt_sk.decrypt(&c2, &glwe_params, table, &mut decrypt_context);

    assert_eq!(
        mod_t.reduce_mul(input1[0], poly_length as ValueT),
        trace_msg[0]
    );

    assert!(trace_msg[1..].iter().all(|&v| v == 0));

    let mut scalar_residue = glwe_params
        .base_q()
        .wrapping_decompose(poly_length as ValueT, t);

    scalar_residue
        .iter_mut()
        .zip(moduli.iter())
        .for_each(|(s, m)| {
            m.reduce_inv_assign(s);
        });

    c1.mul_scalar_assign(&scalar_residue, poly_length, rns_poly_len, &moduli);

    let mut c2: CrtGlwe<Vec<ValueT>> = CrtGlwe::new(c2.0);

    trace_key.trace_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut trace_context,
    );

    let c2 = c2.into_ntt_form(table);

    let trace_msg = dcrt_sk.decrypt(&c2, &glwe_params, table, &mut decrypt_context);

    assert_eq!(input1[0], trace_msg[0]);
    assert!(trace_msg[1..].iter().all(|&v| v == 0));
}

#[test]
fn test_dcrt_glwe_trace() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 1 << 15;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; _] = [1125899906826241, 1125899906629633];
    let moduli = moduli_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(log_n, &moduli).unwrap();

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        dimension,
        poly_length,
        mod_t,
        mod_gamma,
        &moduli,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let moduli_count = glwe_params.cipher_moduli_count();
    let rns_poly_len = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let trace_key = DcrtGlweTraceKey::new(&glev_params, &dcrt_sk, Arc::new(table), &mut rng);
    let table = trace_key.table();

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c2: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut trace_context =
        DcrtGlweTraceContext::new(dimension, poly_length, rns_poly_len, big_uint_poly_len);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    trace_key.trace_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut trace_context,
    );

    let trace_msg = dcrt_sk.decrypt(&c2, &glwe_params, table, &mut decrypt_context);

    assert_eq!(
        mod_t.reduce_mul(input1[0], poly_length as ValueT),
        trace_msg[0]
    );

    assert!(trace_msg[1..].iter().all(|&v| v == 0));

    c1.mul_scalar_assign(
        &glwe_params
            .base_q()
            .wrapping_decompose(poly_length as ValueT, t)
            .iter()
            .zip(moduli.iter())
            .map(|(&n, m)| m.reduce_inv(n))
            .collect::<Vec<_>>(),
        poly_length,
        rns_poly_len,
        &moduli,
    );

    trace_key.trace_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut trace_context,
    );

    let trace_msg = dcrt_sk.decrypt(&c2, &glwe_params, table, &mut decrypt_context);

    assert_eq!(input1[0], trace_msg[0]);
    assert!(trace_msg[1..].iter().all(|&v| v == 0));
}
