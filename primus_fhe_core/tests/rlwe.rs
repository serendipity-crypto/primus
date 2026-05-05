use primus_fhe_core::{
    NttRlweCiphertext, NttRlweSecretKey, RingSecretKeyType, RlweParameters, RlweSecretKey,
};
use primus_integer::UnsignedInteger;
use primus_modulus::BarrettModulus;
use primus_ntt::{NttTable, UintNttTable};
use primus_poly::Polynomial;
use primus_reduce::FieldContext;

const POLY_LENGTH: usize = 1024;
const NOISE_ALPHA: f64 = 2.980_232_238_769_5312e-8; // 2^-25
const SECRET_KEY_GAUSSIAN_STANDARD_DEVIATION: f64 = 3.2;
const PLAIN_MODULI: [usize; 2] = [256, 257];
const SECRET_KEY_TYPES: [RingSecretKeyType; 3] = [
    RingSecretKeyType::Binary,
    RingSecretKeyType::Ternary,
    RingSecretKeyType::Gaussian(SECRET_KEY_GAUSSIAN_STANDARD_DEVIATION),
];

fn from_usize<T: UnsignedInteger>(value: usize) -> T {
    T::try_from(value).unwrap()
}

fn noise_standard_deviation<T, M>(cipher_modulus: M) -> f64
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    let q: f64 = cipher_modulus.value_unchecked().as_into();
    (q * NOISE_ALPHA).max(0.7)
}

fn message_polynomial<T: UnsignedInteger>(plain_modulus: usize) -> Polynomial<Vec<T>> {
    Polynomial::new(
        (0..POLY_LENGTH)
            .map(|index| from_usize(index % plain_modulus))
            .collect(),
    )
}

fn assert_rlwe_secret_key_enc_dec<T, M>(cipher_modulus: M)
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    let mut rng = rand::rng();
    let table = UintNttTable::new(POLY_LENGTH.trailing_zeros(), cipher_modulus).unwrap();
    let noise_standard_deviation = noise_standard_deviation(cipher_modulus);

    for secret_key_type in SECRET_KEY_TYPES {
        for plain_modulus_usize in PLAIN_MODULI {
            let params = RlweParameters::new(
                POLY_LENGTH,
                from_usize(plain_modulus_usize),
                cipher_modulus,
                secret_key_type,
                noise_standard_deviation,
            );
            let secret_key = RlweSecretKey::generate(&params, &mut rng);
            let secret_key = NttRlweSecretKey::from_coeff_secret_key(&secret_key, &table);

            let message = message_polynomial(plain_modulus_usize);
            let mut ciphertext: NttRlweCiphertext<Vec<T>> =
                NttRlweCiphertext::zero(POLY_LENGTH * 2);
            secret_key.encrypt_inplace(&message, &mut ciphertext, &params, &table, &mut rng);

            let decrypted = secret_key.decrypt(&ciphertext, &params, &table);
            assert_eq!(decrypted.as_ref(), message.as_ref());

            let ciphertext = secret_key.encrypt_zeros(&params, &table, &mut rng);
            let decrypted = secret_key.decrypt(&ciphertext, &params, &table);
            assert_eq!(decrypted.as_ref(), vec![T::ZERO; POLY_LENGTH]);

            let ciphertext = secret_key.encrypt_multi_zeros(POLY_LENGTH, &params, &table, &mut rng);
            let decrypted: Vec<T> = secret_key.decrypt_multi_messages(&ciphertext, &params, &table);
            assert_eq!(decrypted, vec![T::ZERO; POLY_LENGTH]);
        }
    }
}

#[test]
fn test_rlwe_secret_key_enc_dec_prime_modulus() {
    assert_rlwe_secret_key_enc_dec::<u64, _>(BarrettModulus::new(1_125_899_906_826_241u64));
    assert_rlwe_secret_key_enc_dec::<u32, _>(BarrettModulus::new(132_120_577u32));
}
