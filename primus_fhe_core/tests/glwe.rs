use pretty_assertions::assert_eq;
use primus_fhe_core::{
    CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweDecryptContext,
    DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::glwe::DcrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{CrtPolynomial, Polynomial};
use primus_reduce::FieldContext;

type ValueT = u64;

const DIMENSION: usize = 2;
const POLY_LENGTH: usize = 512;
const NOISE_STANDARD_DEVIATION: f64 = 3.2;
const SECRET_KEY_GAUSSIAN_STANDARD_DEVIATION: f64 = 3.2;
const PLAIN_MODULI: [ValueT; 3] = [256, 257, 12_289];
const GAMMA_MODULUS: ValueT = 2_305_843_009_213_554_689;
const CIPHER_MODULI: [ValueT; 2] = [1_125_899_906_826_241, 1_125_899_906_629_633];
const SECRET_KEY_TYPES: [RingSecretKeyType; 3] = [
    RingSecretKeyType::Binary,
    RingSecretKeyType::Ternary,
    RingSecretKeyType::Gaussian(SECRET_KEY_GAUSSIAN_STANDARD_DEVIATION),
];

fn message_polynomial(plain_modulus: ValueT) -> Polynomial<Vec<ValueT>> {
    Polynomial::new(
        (0..POLY_LENGTH)
            .map(|index| index as ValueT % plain_modulus)
            .collect(),
    )
}

fn decompose_message<M>(
    message: &Polynomial<Vec<ValueT>>,
    params: &CrtGlweParameters<ValueT, M>,
) -> CrtPolynomial<Vec<ValueT>>
where
    M: FieldContext<ValueT>,
{
    let mut decomposed: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(params.rns_poly_len());
    params.base_q().wrapping_decompose_small_polynomial_inplace(
        message,
        &mut decomposed,
        POLY_LENGTH,
        params.plain_modulus_value(),
    );
    decomposed
}

fn assert_dcrt_glwe_secret_key_enc_dec(secret_key_type: RingSecretKeyType, plain_modulus: ValueT) {
    let mod_t = BarrettModulus::new(plain_modulus);
    let mod_gamma = BarrettModulus::new(GAMMA_MODULUS);
    let moduli = CIPHER_MODULI.map(BarrettModulus::new);
    let table = UintCrtNttTable::new(POLY_LENGTH.trailing_zeros(), &moduli).unwrap();
    let mut rng = rand::rng();

    let params = CrtGlweParameters::new(
        DIMENSION,
        POLY_LENGTH,
        mod_t,
        mod_gamma,
        &moduli,
        secret_key_type,
        NOISE_STANDARD_DEVIATION,
    );

    let secret_key = CrtGlweSecretKey::generate(&params, &mut rng);
    let secret_key = DcrtGlweSecretKey::from_coeff_secret_key(&secret_key, &table);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli.len(), POLY_LENGTH);

    let message = message_polynomial(plain_modulus);
    let decomposed_message = decompose_message(&message, &params);
    let mut ciphertext: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(params.rns_glwe_len());

    secret_key.encrypt_inplace(
        &decomposed_message,
        &mut ciphertext,
        &params,
        &table,
        &mut rng,
    );

    let decrypted = secret_key.decrypt(&ciphertext, &params, &table, &mut decrypt_context);
    assert_eq!(decrypted.as_ref(), message.as_ref());

    let mut ciphertext: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(params.rns_glwe_len());
    secret_key.encrypt_zeros_inplace(&mut ciphertext, &params, &table, &mut rng);

    let decrypted = secret_key.decrypt(&ciphertext, &params, &table, &mut decrypt_context);
    assert_eq!(decrypted.as_ref(), vec![0; POLY_LENGTH]);
}

#[test]
fn test_dcrt_glwe_secret_key_enc_dec_crt_modulus() {
    for secret_key_type in SECRET_KEY_TYPES {
        for plain_modulus in PLAIN_MODULI {
            assert_dcrt_glwe_secret_key_enc_dec(secret_key_type, plain_modulus);
        }
    }
}

#[test]
fn test_dcrt_glwe_secret_key_ciphertext_ops_crt_modulus() {
    let plain_modulus = 12_289;
    let mod_t = BarrettModulus::new(plain_modulus);
    let mod_gamma = BarrettModulus::new(GAMMA_MODULUS);
    let moduli = CIPHER_MODULI.map(BarrettModulus::new);
    let table = UintCrtNttTable::new(POLY_LENGTH.trailing_zeros(), &moduli).unwrap();
    let mut rng = rand::rng();

    let params = CrtGlweParameters::new(
        DIMENSION,
        POLY_LENGTH,
        mod_t,
        mod_gamma,
        &moduli,
        RingSecretKeyType::Ternary,
        NOISE_STANDARD_DEVIATION,
    );

    let rns_poly_len = params.rns_poly_len();
    let rns_glwe_len = params.rns_glwe_len();
    let secret_key = CrtGlweSecretKey::generate(&params, &mut rng);
    let secret_key = DcrtGlweSecretKey::from_coeff_secret_key(&secret_key, &table);
    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli.len(), POLY_LENGTH);

    let m0 = message_polynomial(plain_modulus);
    let mut m1 = Polynomial::new(
        (0..POLY_LENGTH)
            .map(|index| (3 * index as ValueT + 1) % plain_modulus)
            .collect(),
    );
    let m2 = Polynomial::new(
        (0..POLY_LENGTH)
            .map(|index| (index as ValueT + 1) % 2)
            .collect(),
    );

    let msg0 = decompose_message(&m0, &params);
    let msg1 = decompose_message(&m1, &params);
    let msg2 = decompose_message(&m2, &params);

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
    let mut c2: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);

    secret_key.encrypt_inplace(&msg0, &mut c0, &params, &table, &mut rng);
    let mut decrypted = secret_key.decrypt(&c0, &params, &table, &mut decrypt_context);
    assert_eq!(decrypted.as_ref(), m0.as_ref());

    secret_key.encrypt_inplace(&msg1, &mut c1, &params, &table, &mut rng);
    c1.add_element_wise_assign(&c0, POLY_LENGTH, rns_poly_len, &moduli);
    m1.add_assign(&m0, mod_t);

    secret_key.decrypt_inplace(&c1, &mut decrypted, &params, &table, &mut decrypt_context);
    assert_eq!(m1, decrypted);

    c1.sub_element_wise_assign(&c0, POLY_LENGTH, rns_poly_len, &moduli);
    m1.sub_assign(&m0, mod_t);

    secret_key.decrypt_inplace(&c1, &mut decrypted, &params, &table, &mut decrypt_context);
    assert_eq!(m1, decrypted);

    let msg2 = table.transform_inplace(msg2);
    let mut expected_product: Polynomial<Vec<ValueT>> = Polynomial::zero(POLY_LENGTH);

    c1.mul_dcrt_polynomial_inplace(&msg2, &mut c2, POLY_LENGTH, &moduli);
    m1.naive_mul_inplace(&m2, &mut expected_product, mod_t);

    secret_key.decrypt_inplace(&c2, &mut decrypted, &params, &table, &mut decrypt_context);
    assert_eq!(expected_product, decrypted);

    c1.neg_assign(rns_poly_len, POLY_LENGTH, &moduli);
    m1.neg_assign(mod_t);

    secret_key.decrypt_inplace(&c1, &mut decrypted, &params, &table, &mut decrypt_context);
    assert_eq!(m1, decrypted);
}
