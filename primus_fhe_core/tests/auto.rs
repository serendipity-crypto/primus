use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweAutoContext, CrtGlweAutoKey, CrtGlweParameters, CrtGlweSecretKey,
    DcrtGlweCiphertext, DcrtGlweDecryptContext, DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::glwe::{CrtGlwe, DcrtGlwe};
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{Polynomial, crt::CrtPolynomial};

#[test]
fn test_rns_glwe_auto() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(log_n, &qi).unwrap();

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        dimension,
        poly_length,
        mod_t,
        mod_gamma,
        &qi,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let moduli_count = qi.len();
    let big_uint_value_len = glwe_params.big_uint_value_len();
    let crt_poly_length = moduli_count * poly_length;
    let big_uint_poly_len = big_uint_value_len * poly_length;

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    assert_eq!(crt_glwe_len, (dimension + 1) * crt_poly_length);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let crt_glev_len = basis.decompose_length() * crt_glwe_len;
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let auto_key = CrtGlweAutoKey::new(&glev_params, 1, &sk, &dcrt_sk, Arc::new(table), &mut rng);

    let table = auto_key.table();

    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    let input: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

    let mut crt_poly: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
    glwe_params
        .base_q()
        .decompose_small_polynomial_inplace(&input, &mut crt_poly, poly_length);

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);
    dcrt_sk.encrypt_inplace(&crt_poly, &mut c0, &glwe_params, table, &mut rng);

    let c1 = c0.into_coeff_form(table);

    let mut c2: CrtGlwe<Vec<ValueT>> = CrtGlwe::zero(crt_glwe_len);
    let mut context = CrtGlweAutoContext::new(poly_length, crt_poly_length, big_uint_poly_len);
    auto_key.automorphism_inplace(
        &c1,
        &mut c2,
        &glev_params,
        glwe_params.base_q(),
        &mut context,
    );

    let c3 = c2.into_ntt_form(table);

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);
    let msg = dcrt_sk.decrypt(&c3, &glwe_params, table, &mut decrypt_context);

    println!("{:?}", input.as_ref());
    println!("{:?}", msg.as_ref());
}
