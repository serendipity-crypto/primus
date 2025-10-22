use itertools::izip;
use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{BigIntegerOps, DivRemScalar, multiply_many_values};
use primus_modulo::{InvModulo, Modulo};
use primus_modulus::BarrettModulus;
use primus_poly::crt::CrtPolynomial;
use primus_reduce::ops::*;
use primus_rns::{BaseConverter, RNSBase};
use rand::distr::{Distribution, Uniform};

type ValueT = u64;

#[test]
fn test_rns() {
    let moduli = [3, 5, 7].map(BarrettModulus::<ValueT>::new);
    let base = RNSBase::new(&moduli).unwrap();

    let residues = &[2, 3, 2];
    let value = base.compose(residues);
    let dec = base.decompose(&value);
    assert_eq!(dec, residues);

    println!("Result: {:?}", value);
}

#[test]
fn test_rns2() {
    let moduli = [256, 257].map(BarrettModulus::<ValueT>::new);
    let base = RNSBase::new(&moduli).unwrap();

    let residues = &[2, 3];
    let value = base.compose(residues);
    let dec = base.decompose(&value);
    assert_eq!(dec, residues);

    println!("Result: {:?}", value);
}

#[test]
fn test_bfv_dec() {
    let mut rng = rand::rng();

    let moduli_value: [ValueT; 2] = [1099511592961, 1099511590913];
    let moduli_count = moduli_value.len();
    let moduli = moduli_value.map(<BarrettModulus<ValueT>>::new);
    let base_q = RNSBase::new(&moduli).unwrap();
    let q = base_q.moduli_product().to_vec();

    let t = 257;
    let gamma = 1125899906629633;
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
    let inv_gamma_mod_t = ShoupFactor::new(mod_t.reduce_inv(mod_t.reduce(gamma)), t);
    let t_gamma_value = multiply_many_values(&[t, gamma]);
    let t_gamma_mod_q = base_q.decompose(&t_gamma_value);
    let plain_uniform = Uniform::new(0, t).unwrap();

    let input: ValueT = plain_uniform.sample(&mut rng);

    println!("input: {input}");

    let mut delta = vec![0; moduli_count];
    let rem = DivRemScalar::div_rem_scalar(&q, t, &mut delta);
    println!("rem: {rem}");
    let delta_residue = base_q.decompose(&delta);

    let converter = BaseConverter::new(&base_q, &base_t_gamma);

    let mut msg_residue: Vec<ValueT> = vec![0; moduli_count];
    msg_residue
        .iter_mut()
        .zip(&moduli)
        .for_each(|(r, modulus)| {
            *r = modulus.reduce(input);
        });

    // delta * m
    izip!(msg_residue.iter_mut(), delta_residue.iter(), &moduli).for_each(|(a, &b, m)| {
        m.reduce_mul_assign(a, b);
    });

    let mut temp = vec![0; moduli_count];
    let _ = delta.slice_mul_value_inplace(input, &mut temp);
    println!("delta * m: {:?}", base_q.compose(&msg_residue));
    println!("delta * m: {:?}", temp);

    // Add noise
    moduli[0].reduce_add_assign(&mut msg_residue[0], 100);
    moduli[1].reduce_add_assign(&mut msg_residue[1], 100);

    // gamma * t * delta * m
    izip!(msg_residue.iter_mut(), t_gamma_mod_q.iter(), &moduli).for_each(|(a, &b, m)| {
        m.reduce_mul_assign(a, b);
    });

    let mut msg_t_gamma: Vec<ValueT> = vec![0; 2];

    println!("gamma * t * delta * m: {:?}", base_q.compose(&msg_residue));

    converter.fast_convert(msg_residue.as_ref(), msg_t_gamma.as_mut());

    println!(
        "gamma * t * delta * m: {:?}",
        base_t_gamma.compose(&msg_t_gamma)
    );

    izip!(
        msg_t_gamma.iter_mut(),
        minus_inv_q_mod_t_gamma.iter(),
        &t_gamma
    )
    .for_each(|(a, &b, m)| {
        m.reduce_mul_assign(a, b);
    });

    println!("gamma * m: {}", gamma * input);
    let y = base_t_gamma.compose(&msg_t_gamma);
    println!("y: {:?}", y);
    println!(
        "gamma * m residue: {:?}",
        base_t_gamma.decompose(&multiply_many_values(&[gamma, input]))
    );
    println!("y residue: {:?}", msg_t_gamma);

    let (y_t, y_gamma) = (msg_t_gamma[0], msg_t_gamma[1]);

    let dec =
        inv_gamma_mod_t.factor_mul_modulo(mod_t.reduce(mod_gamma.reduce_sub(y_t, y_gamma)), t);

    debug_assert_eq!(input, dec);
}

#[test]
fn test_bfv_dec_array() {
    let mut rng = rand::rng();

    let moduli_value: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let moduli_count = moduli_value.len();
    let moduli = moduli_value.map(<BarrettModulus<ValueT>>::new);
    let base_q = RNSBase::new(&moduli).unwrap();
    let q = base_q.moduli_product().to_vec();

    let t = 12289;
    let gamma = 2305843009213554689;
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
    let inv_gamma_mod_t = ShoupFactor::new(gamma.modulo(mod_t).inv_modulo(mod_t), t);
    let t_gamma_value = multiply_many_values(&[t, gamma]);
    let t_gamma_mod_q = base_q.decompose(&t_gamma_value);
    let plain_uniform = Uniform::new(0, t).unwrap();

    let poly_length = 4096;
    let input: Vec<ValueT> = plain_uniform
        .sample_iter(&mut rng)
        .take(poly_length)
        .collect();
    let mut big_uint_values = vec![0; base_q.big_uint_value_len() * poly_length];
    big_uint_values
        .chunks_exact_mut(base_q.big_uint_value_len())
        .zip(input.iter())
        .for_each(|(a, &b)| {
            a[0] = b;
            // if b <= t / 2 {
            //     a[0] = b;
            // } else {
            //     let _ = q.slice_sub_value_inplace(b, a);
            // }
        });

    let mut delta = vec![0; moduli_count];
    let _rem = DivRemScalar::div_rem_scalar(&q, t, &mut delta);
    let delta_mod_q = base_q.decompose(&delta);

    let converter = BaseConverter::new(&base_q, &base_t_gamma);

    let mut msg_q: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(moduli_count * poly_length);
    let mut msg_t_gamma: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(2 * poly_length);

    base_q.decompose_multiple_values_inplace(&big_uint_values, msg_q.as_mut(), poly_length);

    // delta * m
    msg_q.mul_scalar_assign(&delta_mod_q, poly_length, &moduli);

    // gamma * t * delta * m
    msg_q.mul_scalar_assign(&t_gamma_mod_q, poly_length, &moduli);

    converter.fast_convert_array(msg_q.as_ref(), msg_t_gamma.as_mut(), poly_length);

    msg_t_gamma.mul_scalar_assign(&minus_inv_q_mod_t_gamma, poly_length, &t_gamma);

    let mut result = vec![0; poly_length];

    let (y_t_slices, y_gamma_slices) = msg_t_gamma.as_ref().split_at(poly_length);

    izip!(result.iter_mut(), y_t_slices, y_gamma_slices).for_each(|(res, y_t, y_gamma)| {
        *res = inv_gamma_mod_t.factor_mul_modulo(gamma - y_gamma + y_t, t);
    });

    debug_assert_eq!(input, result);
}
