use primus_fhe_core::{
    CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweDecryptContext,
    DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::glwe::DcrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{Polynomial, crt::CrtPolynomial};

type ValueT = u64;

const PLAIN_MODULUS_VALUE: ValueT = 12289;
const N: usize = 512;

#[test]
fn test_rns_glwe() {
    let qi_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let qi = qi_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(N.trailing_zeros(), &qi).unwrap();
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
