use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweAutoContext, CrtGlweAutoKey, CrtGlweParameters, CrtGlweSecretKey,
    DcrtGlweAutoKey, DcrtGlweCiphertext, DcrtGlweDecryptContext, DcrtGlweSecretKey,
    RingSecretKeyType, crt_poly_auto_inplace, dcrt_poly_ntt_auto_inplace,
};
use primus_lattice::glwe::{CrtGlwe, DcrtGlwe};
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{CrtPolynomial, Polynomial};
use rand::RngExt;

#[test]
fn test_crt_glwe_auto() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    // let gamma: ValueT = 2199023190017;
    let gamma: ValueT = 2305843009213554689;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
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
    let crt_poly_length = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let mut auto_degree = rng.random_range(0..poly_length * 2);
    if auto_degree & 1 == 0 {
        auto_degree |= 1;
    }
    println!("degree: {auto_degree}");
    let auto_key = CrtGlweAutoKey::new(
        &glev_params,
        auto_degree,
        &sk,
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = auto_key.table();

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
    let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut auto_c1: CrtGlwe<Vec<ValueT>> = CrtGlwe::zero(rns_glwe_len);
    let mut c3: CrtGlwe<Vec<ValueT>> = CrtGlwe::zero(rns_glwe_len);
    let mut auto_context = CrtGlweAutoContext::new(poly_length, crt_poly_length, big_uint_poly_len);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    let c1 = c1.into_coeff_form(table);

    c1.iter_crt_poly(crt_poly_length)
        .zip(auto_c1.iter_crt_poly_mut(crt_poly_length))
        .for_each(|(in_crt_poly, auto_crt_poly)| {
            crt_poly_auto_inplace(
                in_crt_poly.0,
                auto_crt_poly.0,
                auto_key.auto_helper(),
                poly_length,
                &moduli,
            );
        });

    let mut auto_sk = sk.clone();
    sk.iter_crt_poly()
        .zip(auto_sk.key_mut().chunks_exact_mut(crt_poly_length))
        .for_each(|(in_crt_poly, auto_crt_poly)| {
            crt_poly_auto_inplace(
                in_crt_poly.0,
                auto_crt_poly,
                auto_key.auto_helper(),
                poly_length,
                &moduli,
            );
        });
    let dcrt_auto_sk = DcrtGlweSecretKey::from_coeff_secret_key(&auto_sk, table);

    let auto_c1 = auto_c1.into_ntt_form(table);
    let auto_msg_1 = dcrt_auto_sk.decrypt(&auto_c1, &glwe_params, table, &mut decrypt_context);

    auto_key.automorphism_inplace(
        &c1,
        &mut c3,
        &glev_params,
        glwe_params.base_q(),
        &mut auto_context,
    );

    let c3 = c3.into_ntt_form(table);

    let auto_msg_2 = dcrt_sk.decrypt(&c3, &glwe_params, table, &mut decrypt_context);

    assert_eq!(auto_msg_1.as_ref(), auto_msg_2.as_ref());
}

#[test]
fn test_dcrt_glwe_auto() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2305843009213554689;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
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
    let crt_poly_length = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let mut auto_degree = rng.random_range(0..poly_length * 2);
    if auto_degree & 1 == 0 {
        auto_degree |= 1;
    }
    println!("degree: {auto_degree}");

    let auto_key = DcrtGlweAutoKey::new(
        &glev_params,
        auto_degree,
        &dcrt_sk,
        Arc::new(table),
        &mut rng,
    );
    let table = auto_key.table();

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
    let mut c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut auto_c1: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c3: DcrtGlweCiphertext<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut auto_context = CrtGlweAutoContext::new(poly_length, crt_poly_length, big_uint_poly_len);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, table, &mut rng);

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, table, &mut decrypt_context);
    assert_eq!(m_dec, input1);

    // Manually apply NTT automorphism to each polynomial of ciphertext
    c1.iter_dcrt_poly(crt_poly_length)
        .zip(auto_c1.iter_dcrt_poly_mut(crt_poly_length))
        .for_each(|(in_dcrt_poly, auto_dcrt_poly)| {
            dcrt_poly_ntt_auto_inplace(
                in_dcrt_poly.0,
                auto_dcrt_poly.0,
                auto_key.auto_helper(),
                poly_length,
            );
        });

    let mut dcrt_auto_sk = dcrt_sk.clone();
    dcrt_sk
        .iter_dcrt_poly()
        .zip(dcrt_auto_sk.iter_dcrt_poly_mut())
        .for_each(|(in_dcrt_poly, auto_dcrt_poly)| {
            dcrt_poly_ntt_auto_inplace(
                in_dcrt_poly.0,
                auto_dcrt_poly.0,
                auto_key.auto_helper(),
                poly_length,
            );
        });

    let auto_msg_1 = dcrt_auto_sk.decrypt(&auto_c1, &glwe_params, table, &mut decrypt_context);

    auto_key.automorphism_inplace(
        &c1,
        &mut c3,
        &glev_params,
        glwe_params.base_q(),
        &mut auto_context,
    );

    let auto_msg_2 = dcrt_sk.decrypt(&c3, &glwe_params, table, &mut decrypt_context);

    assert_eq!(auto_msg_1.as_ref(), auto_msg_2.as_ref());
}
