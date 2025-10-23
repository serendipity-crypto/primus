use indicatif::{ProgressIterator, ProgressStyle};
use pretty_assertions::assert_eq;
use primus_fhe_core::{
    CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweDecryptContext,
    DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::glwe::DcrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{Dcrt, DcrtTable, UintCrtNttTable};
use primus_poly::{Polynomial, crt::CrtPolynomial};

#[test]
fn test_rns_glwe() {
    type ValueT = u64;

    let dimension = 2;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2305843009213554689;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let moduli_count = qi.len();
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

    let crt_poly_length = moduli_count * poly_length;

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    assert_eq!(crt_glwe_len, (dimension + 1) * moduli_count * poly_length);

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    let style = ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos:>7}/{len:7} ({eta})",
    )
    .unwrap()
    .progress_chars("##-");

    for _ in (0..1000).progress_with_style(style) {
        let m0: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

        let mut msg0: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
        glwe_params
            .base_q()
            .decompose_small_polynomial_inplace(&m0, &mut msg0, poly_length);

        let m1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

        let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
        glwe_params
            .base_q()
            .decompose_small_polynomial_inplace(&m1, &mut msg1, poly_length);

        let m2: Polynomial<Vec<ValueT>> = Polynomial::random_binary(poly_length, &mut rng);
        let mut msg2: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(crt_poly_length);
        glwe_params
            .base_q()
            .decompose_small_polynomial_inplace(&m2, &mut msg2, poly_length);

        let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);
        let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);
        let mut c2: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);

        // Encryption and Decryption
        dcrt_sk.encrypt_inplace(&msg0, &mut c0, &glwe_params, &table, &mut rng);

        let mut m_dec = dcrt_sk.decrypt(&c0, &glwe_params, &table, &mut decrypt_context);

        debug_assert_eq!(m0.as_ref(), m_dec.as_ref());

        // Ciphertext Addition
        dcrt_sk.encrypt_inplace(&msg1, &mut c1, &glwe_params, &table, &mut rng);

        c1.add_element_wise_assign(&c0, poly_length, crt_poly_length, &qi);

        dcrt_sk.decrypt_inplace(&c1, &mut m_dec, &glwe_params, &table, &mut decrypt_context);

        assert_eq!(m0.clone().add(&m1, mod_t), m_dec);

        // Ciphertext Subtraction
        c1.sub_element_wise_assign(&c0, poly_length, crt_poly_length, &qi);

        dcrt_sk.decrypt_inplace(&c1, &mut m_dec, &glwe_params, &table, &mut decrypt_context);

        assert_eq!(m1, m_dec);

        // Ciphertext-to-Plaintext Multiplication
        let msg2 = table.transform_inplace(msg2);

        c1.mul_dcrt_polynomial_inplace(&msg2, &mut c2, poly_length, &qi);

        dcrt_sk.decrypt_inplace(&c2, &mut m_dec, &glwe_params, &table, &mut decrypt_context);

        let mut m1_mul_m2: Polynomial<Vec<ValueT>> = Polynomial::zero(poly_length);
        m1.naive_mul_inplace(&m2, &mut m1_mul_m2, mod_t);

        assert_eq!(m1_mul_m2, m_dec);
    }
}
