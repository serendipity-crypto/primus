use primus_fhe_core::{
    CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_integer::{DivRemScalar, izip, multiply_many_values};
use primus_lattice::glwe::DcrtGlwe;
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{BigUintPolynomial, crt::CrtPolynomial};
use primus_reduce::ops::*;
use primus_rns::{BaseConverter, RNSBase};
use rand::distr::{Distribution, Uniform};

type ValueT = u64;

const PLAIN_MODULUS_VALUE: ValueT = 256;
const N: usize = 1024;

#[test]
fn test_rns_glwe() {
    let moduli_value: [ValueT; 2] = [1099511592961, 1099511590913];
    let moduli = moduli_value.map(<BarrettModulus<ValueT>>::new);
    let rns_base = RNSBase::new(&moduli).unwrap();
    let table = UintCrtNttTable::new(N.trailing_zeros(), &moduli).unwrap();
    let modulus = rns_base.moduli_product().to_vec();
    let big_uint_value_len = rns_base.big_uint_value_len();
    let moduli_count = rns_base.moduli_count();

    let plain_uniform = Uniform::new(0, PLAIN_MODULUS_VALUE).unwrap();

    assert_eq!(modulus.len(), 2);

    let mut rng = rand::rng();

    let params = CrtGlweParameters::new(
        2,
        N,
        PLAIN_MODULUS_VALUE,
        &moduli,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let poly_length = params.poly_length();

    let sk = CrtGlweSecretKey::generate(&params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);
    let crt_glwe_len = dcrt_sk.crt_glwe_len();

    let mut c0: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(crt_glwe_len);

    let mut big_uint_poly: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_value_len * poly_length);

    big_uint_poly
        .as_mut_slice()
        .chunks_exact_mut(big_uint_value_len)
        .for_each(|v| v[0] = plain_uniform.sample(&mut rng));

    let mut crt_poly: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);
    let mut msg: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);

    rns_base.decompose_polynomial_inplace(&big_uint_poly, &mut crt_poly, poly_length);

    dcrt_sk.encrypt_inplace(&crt_poly, &mut c0, &params, &table, &mut rng);

    dcrt_sk.phase_inplace(&c0, &mut msg, &params, &table);

    let t = params.plain_modulus_value();
    let gamma = 1125899906826241;

    let t_gamma = [t, gamma].map(<BarrettModulus<ValueT>>::new);
    let base_t_gamma = RNSBase::new(&t_gamma).unwrap();
    let q_mod_t_gamma = base_t_gamma.decompose(params.cipher_modulus());
    let minus_inv_q_mod_t_gamma: Vec<ValueT> = q_mod_t_gamma
        .iter()
        .zip(&t_gamma)
        .map(|(&x, m)| m.reduce_neg(m.reduce_inv(x)))
        .collect();
    let inv_gamma_mod_t = t_gamma[0].reduce(gamma);

    let converter = BaseConverter::new(&rns_base, &base_t_gamma);
    let mut crt_poly_out: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(2 * poly_length);

    converter.fast_convert_array(msg.as_ref(), crt_poly_out.as_mut(), poly_length);

    crt_poly_out.mul_scalar_assign(&minus_inv_q_mod_t_gamma, poly_length, &t_gamma);

    let (y_t, y_gamma) = crt_poly_out.as_ref().split_at(poly_length);
    let dec: Vec<_> = izip!(y_gamma, y_t)
        .map(|(y_t_i, y_gamma_i)| {
            if y_t_i >= y_gamma_i {
                t_gamma[0].reduce_mul(y_t_i - y_gamma_i, inv_gamma_mod_t)
            } else {
                let temp = t_gamma[0].reduce(y_gamma_i - y_t_i);
                t_gamma[0].reduce_mul(t_gamma[0].reduce_neg(temp), inv_gamma_mod_t)
            }
        })
        .collect();

    // rns_base.compose_polynomial_inplace(&msg, &mut big_uint_poly2, poly_length);

    let from: Vec<_> = big_uint_poly
        .as_slice()
        .iter()
        .step_by(big_uint_value_len)
        .copied()
        .collect();

    debug_assert_eq!(from.as_slice(), dec.as_slice());
}

#[test]
fn test_rns() {
    let mut rng = rand::rng();

    let moduli_value: [ValueT; 2] = [1099511592961, 1099511590913];
    let moduli_count = moduli_value.len();
    let moduli = moduli_value.map(<BarrettModulus<ValueT>>::new);
    let base_q = RNSBase::new(&moduli).unwrap();
    let q = base_q.moduli_product().to_vec();
    let poly_length = 8;

    let t = 257;
    let gamma = 1125899906826241;
    let mod_t = <BarrettModulus<ValueT>>::new(t);
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);
    let t_gamma = [mod_t, mod_gamma];
    let base_t_gamma = RNSBase::new(&t_gamma).unwrap();
    let q_mod_t_gamma = base_t_gamma.decompose(&q);
    let minus_inv_q_mod_t_gamma: Vec<ValueT> = q_mod_t_gamma
        .iter()
        .zip(&t_gamma)
        .map(|(&x, m)| m.reduce_neg(m.reduce_inv(x)))
        .collect();
    let inv_gamma_mod_t = mod_t.reduce(gamma);
    let t_gamma_value = multiply_many_values(&[t, gamma]);
    let t_gamma_mod_q = base_q.decompose(&t_gamma_value);
    let plain_uniform = Uniform::new(0, t).unwrap();

    let input: Vec<ValueT> = plain_uniform
        .sample_iter(&mut rng)
        .take(poly_length)
        .collect();

    let mut delta = vec![0, 0];
    let _rem = DivRemScalar::div_rem_scalar(&q, t, &mut delta);
    let delta_residue = base_q.decompose(&delta);

    let converter = BaseConverter::new(&base_q, &base_t_gamma);

    let mut msg: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);
    msg.iter_each_modulus_mut(poly_length)
        .zip(&moduli)
        .for_each(|(m, modulus)| {
            m.iter_mut().zip(input.iter()).for_each(|(a, b)| {
                *a = modulus.reduce(*b);
            });
        });
    msg.mul_scalar_assign(&delta_residue, poly_length, &moduli);

    msg.mul_scalar_assign(&t_gamma_mod_q, poly_length, &moduli);

    let mut crt_poly_out: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(2 * poly_length);

    converter.fast_convert_array(msg.as_ref(), crt_poly_out.as_mut(), poly_length);

    crt_poly_out.mul_scalar_assign(&minus_inv_q_mod_t_gamma, poly_length, &t_gamma);

    let (y_t, y_gamma) = crt_poly_out.as_ref().split_at(poly_length);
    let dec: Vec<_> = izip!(y_gamma, y_t)
        .map(|(y_t_i, y_gamma_i)| {
            if y_t_i >= y_gamma_i {
                t_gamma[0].reduce_mul(y_t_i - y_gamma_i, inv_gamma_mod_t)
            } else {
                let temp = t_gamma[0].reduce(y_gamma_i - y_t_i);
                t_gamma[0].reduce_mul(t_gamma[0].reduce_neg(temp), inv_gamma_mod_t)
            }
        })
        .collect();

    debug_assert_eq!(input.as_slice(), dec.as_slice());
}
