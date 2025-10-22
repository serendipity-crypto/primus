use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweAutoContext, CrtGlweAutoKey, CrtGlweParameters, CrtGlweSecretKey,
    DcrtGlweCiphertext, DcrtGlweDecryptContext, DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_integer::BigIntegerOps;
use primus_lattice::{
    context::DcrtGlevContext,
    glev::DcrtGlev,
    glwe::{CrtGlwe, DcrtGlwe},
};
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, Ntt, NttTable, UintCrtNttTable, UintNttTable};
use primus_poly::{BigUintPolynomial, Polynomial, crt::CrtPolynomial};

type ValueT = u64;

const PLAIN_MODULUS_VALUE: ValueT = 12289;
const N: usize = 512;
const LOG_N: u32 = N.trailing_zeros();

#[test]
fn test_rns_glwe() {
    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(LOG_N, &qi).unwrap();
    let mod_t = <BarrettModulus<ValueT>>::new(PLAIN_MODULUS_VALUE);

    let mut rng = rand::rng();

    let params = CrtGlweParameters::new(
        2,
        N,
        mod_t,
        <BarrettModulus<ValueT>>::new(2305843009213554689),
        &qi,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let poly_length = params.poly_length();

    let sk = CrtGlweSecretKey::generate(&params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    let input: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

    let mut crt_poly: CrtPolynomial<Vec<ValueT>> =
        CrtPolynomial::zero(params.cipher_moduli_count() * poly_length);
    params
        .base_q()
        .decompose_small_polynomial_inplace(&input, &mut crt_poly, poly_length);

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);
    dcrt_sk.encrypt_inplace(&crt_poly, &mut c0, &params, &table, &mut rng);

    let mut decrypt_context =
        DcrtGlweDecryptContext::new(params.cipher_moduli_count(), poly_length);
    let msg = dcrt_sk.decrypt(&c0, &params, &table, &mut decrypt_context);

    debug_assert_eq!(input.as_ref(), msg.as_ref());
}

#[test]
fn test_poly_mul_rns_glev() {
    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(LOG_N, &qi).unwrap();
    let mod_t = <BarrettModulus<ValueT>>::new(PLAIN_MODULUS_VALUE);

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        2,
        N,
        mod_t,
        <BarrettModulus<ValueT>>::new(2305843009213554689),
        &qi,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let dimension = glwe_params.dimension();
    let poly_length = glwe_params.poly_length();
    let moduli_count = glwe_params.cipher_moduli_count();
    let big_uint_value_len = glwe_params.big_uint_value_len();
    let crt_poly_len = poly_length * moduli_count;
    let big_uint_poly_len = poly_length * big_uint_value_len;

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = (dimension + 1) * crt_poly_len;

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 2, None, glwe_params.base_q());
    let decompose_length = basis.decompose_length();
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);
    let crt_glev_len = decompose_length * crt_glwe_len;

    let mut dcrt_glev: DcrtGlev<Vec<ValueT>> = DcrtGlev::zero(crt_glev_len);

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_len);
    glwe_params
        .base_q()
        .decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length);

    dcrt_sk.encrypt_dcrt_glev_inplace(&msg1, &mut dcrt_glev, &glev_params, &table, &mut rng);

    let input2: Polynomial<Vec<ValueT>> = Polynomial::random_binary(poly_length, &mut rng);
    let mut msg2_big_poly: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_poly_len);
    msg2_big_poly
        .iter_mut(big_uint_value_len)
        .zip(input2.iter())
        .for_each(|(a, &b)| {
            if b <= PLAIN_MODULUS_VALUE / 2 {
                a[0] = b;
            } else {
                let _ = glwe_params.cipher_modulus().slice_sub_value_inplace(b, a);
            }
        });

    let mut result: DcrtGlwe<Vec<ValueT>> = DcrtGlwe::zero(crt_glwe_len);

    let mut context = DcrtGlevContext::new(poly_length, crt_poly_len, big_uint_poly_len);

    dcrt_glev.mul_big_uint_poly_inplace(
        &msg2_big_poly,
        &mut result,
        glev_params.basis(),
        &table,
        glwe_params.base_q(),
        &mut context,
    );

    let mut context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    let dec = dcrt_sk.decrypt(&result, &glwe_params, &table, &mut context);

    let ntt_table = UintNttTable::new(LOG_N, mod_t).unwrap();
    let input1_ntt = ntt_table.transform_inplace(input1.clone());
    let input2_ntt = ntt_table.transform_inplace(input2.clone());
    let dec_ntt = ntt_table.transform_inplace(dec.clone());

    let mul = input1_ntt.mul(&input2_ntt, mod_t);
    assert_eq!(dec_ntt, mul);
}

#[test]
fn test_rns_glwe_auto() {
    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let table = Arc::new(UintCrtNttTable::new(LOG_N, &qi).unwrap());
    let mod_t = <BarrettModulus<ValueT>>::new(PLAIN_MODULUS_VALUE);

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        2,
        N,
        mod_t,
        <BarrettModulus<ValueT>>::new(2305843009213554689),
        &qi,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let poly_length = glwe_params.poly_length();
    let moduli_count = glwe_params.cipher_moduli_count();
    let crt_poly_len = poly_length * moduli_count;
    let big_uint_poly_len = poly_length * glwe_params.big_uint_value_len();

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 2, None, glwe_params.base_q());
    let decompose_length = basis.decompose_length();
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, table.as_ref());

    let auto_key = CrtGlweAutoKey::new(
        &glev_params,
        poly_length + 1,
        &sk,
        &dcrt_sk,
        Arc::clone(&table),
        &mut rng,
    );

    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    let input: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

    let mut crt_poly: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_len);
    glwe_params
        .base_q()
        .decompose_small_polynomial_inplace(&input, &mut crt_poly, poly_length);

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);
    dcrt_sk.encrypt_inplace(&crt_poly, &mut c0, &glwe_params, table.as_ref(), &mut rng);

    let c1 = c0.into_coeff_form(table.as_ref());

    let mut c2: CrtGlwe<Vec<ValueT>> = CrtGlwe::zero(crt_glwe_len);
    let mut context = CrtGlweAutoContext::new(poly_length, crt_poly_len, big_uint_poly_len);
    auto_key.automorphism_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut context,
    );

    let c3 = c2.into_ntt_form(table.as_ref());

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);
    let msg = dcrt_sk.decrypt(&c3, &glwe_params, table.as_ref(), &mut decrypt_context);

    println!("{:?}", input.as_ref());
    println!("{:?}", msg.as_ref());
}
