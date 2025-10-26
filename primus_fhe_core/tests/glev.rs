use itertools::izip;
use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweParameters, CrtGlweSecretKey, DcrtGlweDecryptContext,
    DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::{context::DcrtGlevContext, glev::DcrtGlev, glwe::DcrtGlwe};
use primus_modulus::BarrettModulus;
use primus_ntt::{Dcrt, DcrtTable, UintCrtNttTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Polynomial, crt::CrtPolynomial, dcrt::DcrtPolynomial,
};

#[test]
fn test_rns_glev() {
    type ValueT = u64;

    let dimension = 3;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; 2] = [1125899906826241, 1125899906629633];
    let moduli = moduli_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(log_n, &moduli).unwrap();

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        dimension,
        poly_length,
        mod_t,
        mod_gamma,
        &moduli,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let rns_glwe_len = glwe_params.rns_glwe_len();
    let moduli_count = glwe_params.cipher_moduli_count();
    let rns_poly_len = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);
    let rns_glev_len = glev_params.rns_glev_len();

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

    let mut glev_context = DcrtGlevContext::new(poly_length, rns_poly_len, big_uint_poly_len);

    let mut dcrt_glev: DcrtGlev<Vec<ValueT>> = DcrtGlev::zero(rns_glev_len);

    let mut desired: Polynomial<Vec<ValueT>> = Polynomial::zero(poly_length);

    let input1: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
    let input2: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);

    input1.naive_mul_inplace(&input2, &mut desired, mod_t);

    let mut msg2_big_uint_poly: BigUintPolynomial<Vec<ValueT>> =
        BigUintPolynomial::zero(big_uint_poly_len);

    let mut msg1: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
    let mut msg2: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input1, &mut msg1, poly_length, t);

    glwe_params
        .base_q()
        .wrapping_decompose_small_polynomial_inplace(&input2, &mut msg2, poly_length, t);

    msg2.mul_scalar_assign(glwe_params.delta_mod_q(), poly_length, &moduli);

    glwe_params
        .base_q()
        .compose_polynomial_inplace(&msg2, &mut msg2_big_uint_poly, poly_length);

    let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlwe::zero(rns_glwe_len);

    dcrt_sk.encrypt_dcrt_glev_inplace(&msg1, &mut dcrt_glev, &glev_params, &table, &mut rng);

    dcrt_glev.mul_big_uint_poly_inplace(
        &msg2_big_uint_poly,
        &mut c1,
        glev_params.basis(),
        &table,
        glwe_params.base_q(),
        &mut glev_context,
    );

    // c1.add_dcrt_glev_mul_big_uint_poly_assign(
    //     &dcrt_glev,
    //     &msg2_big_uint_poly,
    //     glev_params.basis(),
    //     &table,
    //     glwe_params.base_q(),
    //     &mut glev_context,
    // );

    let m_dec = dcrt_sk.decrypt(&c1, &glwe_params, &table, &mut decrypt_context);

    pretty_assertions::assert_eq!(m_dec, desired);
}

#[test]
fn test_key_switching() {
    type ValueT = u64;

    let dimension = 3;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    let mod_gamma = <BarrettModulus<ValueT>>::new(gamma);

    let moduli_values: [ValueT; _] = [1125899906826241, 1125899906629633];
    let moduli = moduli_values.map(<BarrettModulus<ValueT>>::new);
    let table = UintCrtNttTable::new(log_n, &moduli).unwrap();

    let mut rng = rand::rng();

    let glwe_params = CrtGlweParameters::new(
        dimension,
        poly_length,
        mod_t,
        mod_gamma,
        &moduli,
        RingSecretKeyType::Ternary,
        3.20,
    );

    let basis =
        BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 20, None, glwe_params.base_q());
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let rns_poly_len = glwe_params.rns_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();
    let moduli_count = glwe_params.cipher_moduli_count();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glev_len = glev_params.rns_glev_len();
    let uniform_distrs = glev_params.cipher_moduli_uniform_distr();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let mut dcrt_glevs: Vec<DcrtGlev<Vec<ValueT>>> = (0..dimension)
        .map(|_| DcrtGlev::zero(rns_glev_len))
        .collect();
    let mut msgs: Vec<CrtPolynomial<Vec<ValueT>>> = (0..dimension)
        .map(|_| CrtPolynomial::zero(rns_poly_len))
        .collect();

    sk.iter_crt_poly()
        .zip(msgs.iter_mut())
        .for_each(|(a, b)| b.as_mut().copy_from_slice(a));

    msgs.iter()
        .zip(dcrt_glevs.iter_mut())
        .for_each(|(msg, glev)| {
            dcrt_sk.encrypt_dcrt_glev_inplace(&msg, glev, &glev_params, &table, &mut rng);
        });

    let mut cipher: Vec<DcrtPolynomial<Vec<ValueT>>> = (0..dimension)
        .map(|_| DcrtPolynomial::zero(rns_poly_len))
        .collect();

    let mut b: DcrtPolynomial<Vec<ValueT>> = DcrtPolynomial::zero(rns_poly_len);

    primus_distr::sample_crt_gaussian_values_inplace(
        b.as_mut(),
        poly_length,
        &moduli_values,
        glwe_params.noise_distribution(),
        &mut rng,
    );

    table.transform_slice(b.as_mut());

    cipher.iter_mut().for_each(|ai| {
        primus_distr::sample_crt_uniform_values_inplace(
            ai.as_mut(),
            poly_length,
            uniform_distrs,
            &mut rng,
        );
    });

    dcrt_sk
        .iter_dcrt_poly()
        .zip(cipher.iter())
        .for_each(|(si, ai)| {
            b.add_mul_assign(ai, &DcrtPolynomial(ArrayBase(si)), poly_length, &moduli);
        });

    let cipher: Vec<_> = cipher
        .into_iter()
        .map(|a| table.inverse_transform_inplace(a))
        .collect();

    let mut cs: Vec<DcrtGlwe<Vec<ValueT>>> = (0..dimension)
        .map(|_| DcrtGlwe::zero(rns_glwe_len))
        .collect();

    let mut glev_context = DcrtGlevContext::new(poly_length, rns_poly_len, big_uint_poly_len);
    izip!(dcrt_glevs.iter(), cipher.iter(), cs.iter_mut()).for_each(|(glev, ai, result)| {
        glev.mul_crt_poly_inplace(
            ai,
            result,
            glev_params.basis(),
            &table,
            glwe_params.base_q(),
            &mut glev_context,
        );
    });

    let mut res: DcrtGlwe<Vec<ValueT>> = DcrtGlwe::zero(rns_glwe_len);

    let (_, b_) = res.a_b_mut_slices(glev_params.rns_glwe_mid());
    b_.copy_from_slice(b.as_ref());

    let result = cs.iter().fold(res, |mut acc, x| {
        acc.sub_element_wise_assign(x, poly_length, rns_poly_len, &moduli);
        acc
    });

    let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);
    let m_dec = dcrt_sk.decrypt(&result, &glwe_params, &table, &mut decrypt_context);

    println!("{:?}", m_dec.as_ref());
}
