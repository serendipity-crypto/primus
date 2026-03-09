use indicatif::{ProgressIterator, ProgressStyle};
use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_fhe_core::{
    CrtGlevParameters, CrtGlweParameters, CrtGlweSecretKey, DcrtGlweCiphertext,
    DcrtGlweDecryptContext, DcrtGlwePublicKey, DcrtGlweSecretKey, RingSecretKeyType,
};
use primus_lattice::{context::DcrtGlevContext, glwe::DcrtGlwe};
use primus_modulus::BarrettModulus;
use primus_ntt::{DcrtTable, UintCrtNttTable};
use primus_poly::{CrtPolynomial, Polynomial};
use primus_reduce::ops::ReduceNegAssign;
use rand::RngExt;

#[test]
fn test_external_product() {
    type ValueT = u64;

    let dimension = 8;
    let poly_length: usize = 512;
    let log_n = poly_length.trailing_zeros();

    // let t: ValueT = 1 << 15;
    let t: ValueT = 12289;
    let mod_t = <BarrettModulus<ValueT>>::new(t);

    let gamma: ValueT = 2199023190017;
    // let gamma: ValueT = 2305843009213554689;
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

    let moduli_count = glwe_params.cipher_moduli_count();
    let rns_poly_len = glwe_params.rns_poly_len();
    let big_uint_poly_len = glwe_params.big_uint_poly_len();
    let rns_glwe_len = glwe_params.rns_glwe_len();
    let base_q = glwe_params.base_q();

    let sk = CrtGlweSecretKey::generate(&glwe_params, &mut rng);
    let dcrt_sk = DcrtGlweSecretKey::from_coeff_secret_key(&sk, &table);

    let basis = BigUintApproxSignedBasis::new(glwe_params.cipher_modulus(), 30, None, base_q);
    let glev_params = CrtGlevParameters::with_glwe_params(&glwe_params, basis);

    let pk = DcrtGlwePublicKey::new(&dcrt_sk, &glwe_params, &table, &mut rng);

    let style = ProgressStyle::with_template(
        "[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos:>7}/{len:7} ({eta})",
    )
    .unwrap()
    .progress_chars("##-");

    for _ in (0..20).progress_with_style(style) {
        let degree = rng.random_range(0..poly_length);
        // println!("degree: {degree}");

        let ggsw = pk.encrypt_monomial_ggsw(&[1, 1], degree, &glev_params, &table, &mut rng);

        let input: Polynomial<Vec<ValueT>> = Polynomial::random(poly_length, mod_t, &mut rng);
        let mut msg: CrtPolynomial<Vec<ValueT>> = CrtPolynomial::zero(rns_poly_len);
        let mut c1: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);
        let mut c2: DcrtGlwe<Vec<ValueT>> = DcrtGlweCiphertext::zero(rns_glwe_len);

        let mut glev_context = DcrtGlevContext::new(poly_length, rns_poly_len, big_uint_poly_len);
        let mut decrypt_context = DcrtGlweDecryptContext::new(moduli_count, poly_length);

        base_q.wrapping_decompose_small_polynomial_inplace(&input, &mut msg, poly_length, t);

        dcrt_sk.encrypt_inplace(&msg, &mut c1, &glwe_params, &table, &mut rng);

        let c1 = c1.into_coeff_form(&table);

        c1.mul_dcrt_ggsw_inplace(
            &ggsw,
            &mut c2,
            glev_params.basis(),
            &table,
            base_q,
            &mut glev_context,
        );

        let mut input_rt = input.clone();
        input_rt.as_mut_slice().rotate_right(degree);
        input_rt.as_mut_slice()[..degree].iter_mut().for_each(|c| {
            mod_t.reduce_neg_assign(c);
        });

        let output = dcrt_sk.decrypt(&c2, &glwe_params, &table, &mut decrypt_context);

        // input_rt
        //     .iter()
        //     .zip(output.iter())
        //     .enumerate()
        //     .for_each(|(i, (&left, &right))| {
        //         if left != right {
        //             println!("index: {i:4}\nleft: {left:5}\nright:{right:5}\n");
        //         }
        //     });
        assert_eq!(input_rt.as_ref(), output.as_ref());
    }
}
