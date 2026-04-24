use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweExpandCoeffContext, CrtGlweExpandCoeffKey,
    CrtGlweExpandCoeffSyncPool, CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext,
    DcrtGlweDecryptContext, DcrtGlweExpandCoeffContext, DcrtGlweExpandCoeffKey,
    DcrtGlweExpandCoeffSyncPool, DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::glwe::CrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{CrtPolynomial, Polynomial};

#[test]
fn test_crt_glwe_expand_coefficients() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
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

    let expand_key = CrtGlweExpandCoeffKey::new(
        &glev_params,
        glwe_params.base_q(),
        &sk,
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = expand_key.table();

    let mut input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c_expand: Vec<CrtGlwe<Vec<ValueT>>> = vec![CrtGlwe::zero(rns_glwe_len); poly_length];
    let mut expand_context = CrtGlweExpandCoeffContext::new(
        dimension,
        poly_length,
        rns_poly_len,
        big_uint_poly_len,
        moduli_count,
    );
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    let c1 = c1.into_coeff_form(table);

    expand_key.expand_coefficients_inplace(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &mut expand_context,
    );

    for (cipher, &input) in c_expand.into_iter().zip(input1.iter()) {
        let cipher = cipher.into_ntt_form(table);
        let m_dec = dcrt_sk.decrypt(&cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }

    let mut c1 = DcrtGlweCiphertext::new(c1.0);

    input1[256..].fill(0);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let c1 = c1.into_coeff_form(table);

    let mut c_expand: Vec<CrtGlwe<Vec<ValueT>>> = vec![CrtGlwe::zero(rns_glwe_len); 256];

    expand_key.expand_partial_coefficients_inplace(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &mut expand_context,
    );

    for (cipher, &input) in c_expand.into_iter().zip(input1.iter()) {
        let cipher = cipher.into_ntt_form(table);
        let m_dec = dcrt_sk.decrypt(&cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }
}

#[test]
fn test_dcrt_glwe_expand_coefficients() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 12289;
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

    let expand_key = DcrtGlweExpandCoeffKey::new(
        &glev_params,
        glwe_params.base_q(),
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = expand_key.table();

    let mut input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c_expand: Vec<DcrtGlweCiphertext<Vec<ValueT>>> =
        vec![DcrtGlweCiphertext::zero(rns_glwe_len); poly_length];
    let mut expand_context = DcrtGlweExpandCoeffContext::new(
        dimension,
        poly_length,
        rns_poly_len,
        big_uint_poly_len,
        moduli_count,
    );
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    // Input stays in NTT domain for DcrtGlweExpandCoeffKey
    expand_key.expand_coefficients_inplace(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &mut expand_context,
    );

    // Results are already in NTT domain — decrypt directly
    for (cipher, &input) in c_expand.iter().zip(input1.iter()) {
        let m_dec = dcrt_sk.decrypt(cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }

    // Test partial expansion
    input1[256..].fill(0);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&mut msg1, &mut c1, &glwe_params, table, &mut rng);

    let mut c_expand: Vec<DcrtGlweCiphertext<Vec<ValueT>>> =
        vec![DcrtGlweCiphertext::zero(rns_glwe_len); 256];

    expand_key.expand_partial_coefficients_inplace(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &mut expand_context,
    );

    for (cipher, &input) in c_expand.iter().zip(input1.iter()) {
        let m_dec = dcrt_sk.decrypt(cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }
}

#[test]
fn test_dcrt_glwe_expand_coefficients_parallel() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 12289;
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

    let expand_key = DcrtGlweExpandCoeffKey::new(
        &glev_params,
        glwe_params.base_q(),
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = expand_key.table();

    let context_pool = DcrtGlweExpandCoeffSyncPool::with_capacity(
        poly_length.trailing_zeros() as usize,
        dimension,
        poly_length,
        rns_poly_len,
        big_uint_poly_len,
        moduli_count,
    );

    let mut input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c_expand: Vec<DcrtGlweCiphertext<Vec<ValueT>>> =
        vec![DcrtGlweCiphertext::zero(rns_glwe_len); poly_length];
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    // Full parallel expansion
    expand_key.expand_coefficients_inplace_parallel(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &context_pool,
    );

    for (cipher, &input) in c_expand.iter().zip(input1.iter()) {
        let m_dec = dcrt_sk.decrypt(cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }

    // Test partial parallel expansion
    input1[256..].fill(0);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&mut msg1, &mut c1, &glwe_params, table, &mut rng);

    let mut c_expand: Vec<DcrtGlweCiphertext<Vec<ValueT>>> =
        vec![DcrtGlweCiphertext::zero(rns_glwe_len); 256];

    expand_key.expand_partial_coefficients_inplace_parallel(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &context_pool,
    );

    for (cipher, &input) in c_expand.iter().zip(input1.iter()) {
        let m_dec = dcrt_sk.decrypt(cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }
}

#[test]
fn test_crt_glwe_expand_coefficients_parallel() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 12289;
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

    let expand_key = CrtGlweExpandCoeffKey::new(
        &glev_params,
        glwe_params.base_q(),
        &sk,
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = expand_key.table();

    let context_pool = CrtGlweExpandCoeffSyncPool::with_capacity(
        poly_length.trailing_zeros() as usize,
        dimension,
        poly_length,
        rns_poly_len,
        big_uint_poly_len,
        moduli_count,
    );

    let mut input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c_expand: Vec<CrtGlwe<Vec<ValueT>>> = vec![CrtGlwe::zero(rns_glwe_len); poly_length];
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    let c1 = c1.into_coeff_form(table);

    // Full parallel expansion
    expand_key.expand_coefficients_inplace_parallel(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &context_pool,
    );

    for (cipher, &input) in c_expand.into_iter().zip(input1.iter()) {
        let cipher = cipher.into_ntt_form(table);
        let m_dec = dcrt_sk.decrypt(&cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }

    // Partial parallel expansion
    let mut c1 = DcrtGlweCiphertext::new(c1.0);

    input1[256..].fill(0);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let c1 = c1.into_coeff_form(table);

    let mut c_expand: Vec<CrtGlwe<Vec<ValueT>>> = vec![CrtGlwe::zero(rns_glwe_len); 256];

    expand_key.expand_partial_coefficients_inplace_parallel(
        &c1,
        &mut c_expand,
        &glev_params,
        glwe_params.base_q(),
        &context_pool,
    );

    for (cipher, &input) in c_expand.into_iter().zip(input1.iter()) {
        let cipher = cipher.into_ntt_form(table);
        let m_dec = dcrt_sk.decrypt(&cipher, &glwe_params, table, &mut decrypt_context);
        assert_eq!(input, m_dec[0]);
        assert!(m_dec[1..].iter().all(|&v| v == 0));
    }
}
