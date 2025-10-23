use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweParameters, CrtGlweSecretKey, DcrtGlweDecryptContext,
    DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::{context::DcrtGlevContext, glev::DcrtGlev, glwe::DcrtGlwe};
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{BigUintPolynomial, Polynomial, crt::CrtPolynomial};

#[test]
fn test_rns_glev() {
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

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    let mut glev_context = DcrtGlevContext::new(poly_length, crt_poly_length, big_uint_poly_len);

    let mut dcrt_glev: DcrtGlev<Vec<ValueT>> = DcrtGlev::zero(crt_glev_len);

    let mut desired: Polynomial<Vec<ValueT>> = Polynomial::zero(poly_length);

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let input2: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

    input1.naive_mul_inplace(&input2, &mut desired, mod_t);

    let mut msg2_big_uint_poly: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_poly_len);

    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
    let mut msg2: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input2, &mut msg2, poly_length, t);

    msg2.mul_scalar_assign(glwe_params.delta_mod_q(), poly_length, &qi);

    glwe_params
        .base_q()
        .compose_polynomial_inplace(&msg2, &mut msg2_big_uint_poly, poly_length);

    let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlwe::zero(crt_glwe_len);

    dcrt_sk.encrypt_dcrt_glev_inplace(&msg1, &mut dcrt_glev, &glev_params, &table, &mut rng);

    dcrt_glev.mul_big_uint_poly_inplace(
        &msg2_big_uint_poly,
        &mut c1,
        glev_params.basis(),
        &table,
        glwe_params.base_q(),
        &mut glev_context,
    );

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, &table, &mut decrypt_context);

    pretty_assertions::assert_eq!(m_dec, desired);
}
