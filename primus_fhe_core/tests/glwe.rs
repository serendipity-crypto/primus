use primus_fhe_core::{CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweSecretKey};
use primus_lattice::glwe::DcrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{BigUintPolynomial, crt::CrtPolynomial};
use primus_rns::RNSBase;
use rand::distr::{Distribution, Uniform};

type ValueT = u64;
type WideT = u128;

const PLAIN_MODULUS_VALUE: ValueT = 256;
const N: usize = 1024;

#[test]
fn test_crt_glwe_and_dcrt_glwe() {
    let mut rng = rand::rng();
    let moduli_values: [ValueT; 2] = [1099511592961, 1099511590913];
    let moduli = moduli_values.map(<BarrettModulus<ValueT>>::new);
    let rns_base = RNSBase::new(&moduli).unwrap();
    let table = UintCrtNttTable::new(N.trailing_zeros(), &moduli).unwrap();
    let modulus = rns_base.moduli_product().to_vec();
    let big_uint_value_len = modulus.len();
    let moduli_count = moduli.len();

    let plain_uniform = Uniform::new(0, PLAIN_MODULUS_VALUE).unwrap();

    assert_eq!(modulus.len(), 2);

    let modulus_wide = modulus[0] as WideT + ((modulus[1] as WideT) << ValueT::BITS);

    let params = CrtGlweParameters::new(
        2,
        N,
        PLAIN_MODULUS_VALUE,
        &moduli,
        primus_fhe_core::RingSecretKeyType::Ternary,
        3.20,
    );

    let poly_length = params.poly_length();

    let sk = CrtGlweSecretKey::generate(&params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);

    let mut big_uint_poly: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_value_len * poly_length);

    let mut big_uint_poly2: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_value_len * poly_length);

    big_uint_poly
        .as_mut_slice()
        .chunks_exact_mut(big_uint_value_len)
        .for_each(|v| v[0] = plain_uniform.sample(&mut rng));

    let mut crt_poly: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);
    let mut msg: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);

    rns_base.decompose_polynomial_inplace(&big_uint_poly, &mut crt_poly, poly_length);

    dcrt_sk.encrypt_inplace(&crt_poly, &mut c0, &params, &table, &mut rng);

    dcrt_sk.decrypt_inplace(&c0, &mut msg, &params, &table);

    rns_base.compose_polynomial_inplace(&msg, &mut big_uint_poly2, poly_length);

    debug_assert_eq!(big_uint_poly.as_slice(), big_uint_poly2.as_slice());
}
