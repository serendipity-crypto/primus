use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_integer::{UnsignedInteger, izip};
use primus_lattice::{context::DcrtGlevContext, glev::DcrtGlev};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Data, DataMut, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_rns::RNSBase;

use crate::{
    CrtGlevParameters, CrtGlweCiphertext, CrtGlweParameters, DcrtGlweCiphertext, DcrtGlweSecretKey,
};

pub struct CrtGlweKeySwitchingKey<T: UnsignedInteger> {
    key: Vec<T>,
    poly_length: usize,
    crt_poly_length: usize,
    big_uint_poly_len: usize,
    input_crt_glwe_mid: usize,
    dcrt_glev_length: usize,
    output_crt_glwe_mid: usize,
}

impl<T: UnsignedInteger> CrtGlweKeySwitchingKey<T> {
    pub fn new<R, M, Table>(
        input_sk: &DcrtGlweSecretKey<T>,
        input_params: &CrtGlweParameters<T, M>,
        output_sk: &DcrtGlweSecretKey<T>,
        output_params: &CrtGlweParameters<T, M>,
        ksk_params: &CrtGlevParameters<T, M>,
        table: &Table,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        debug_assert_eq!(input_params.poly_length(), output_params.poly_length());
        debug_assert_eq!(
            input_params.cipher_modulus(),
            output_params.cipher_modulus()
        );

        let decompose_length = ksk_params.basis().decompose_length();
        let dcrt_glev_length = decompose_length * output_sk.crt_glwe_len();
        let mut key: Vec<T> = vec![T::ZERO; input_sk.dimension() * dcrt_glev_length];

        izip!(
            key.chunks_exact_mut(dcrt_glev_length),
            input_sk.iter_dcrt_poly(),
        )
        .for_each(|(dcrt_glev, si)| {
            output_sk.encrypt_dcrt_glev_inplace(
                &CrtPolynomial(ArrayBase(si)),
                &mut DcrtGlev::new(ArrayBase(dcrt_glev)),
                &ksk_params,
                table,
                rng,
            );
        });

        let poly_length = input_params.poly_length();
        let crt_poly_length = poly_length * input_params.cipher_moduli_count();
        let big_uint_poly_len = poly_length * input_params.big_uint_value_len();
        let input_crt_glwe_mid = crt_poly_length * input_params.dimension();
        let output_crt_glwe_mid = crt_poly_length * output_params.dimension();
        Self {
            key,
            poly_length,
            crt_poly_length,
            big_uint_poly_len,
            input_crt_glwe_mid,
            dcrt_glev_length,
            output_crt_glwe_mid,
        }
    }

    pub fn key_swithching<M, Table, A, B>(
        &self,
        c_in: &CrtGlweCiphertext<A>,
        c_out: &mut DcrtGlweCiphertext<B>,
        basis: &BigUintApproxSignedBasis<T>,
        table: &Table,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweKeySwitchingContext<T>,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let (a_in, b_in) = c_in.a_b_slices(self.input_crt_glwe_mid);

        let (big_uint_poly, crt_poly, glev_context) = context.as_mut();

        izip!(
            a_in.chunks_exact(self.crt_poly_length),
            self.key.chunks_exact(self.dcrt_glev_length)
        )
        .for_each(|(ai, ki)| {
            let ai = CrtPolynomial(ArrayBase(ai));
            let ki = DcrtGlev::new(ArrayBase(ki));

            rns_base.compose_polynomial_inplace(&ai, big_uint_poly, self.poly_length);

            c_out.add_dcrt_glev_mul_big_uint_poly_assign(
                &ki,
                &big_uint_poly,
                basis,
                table,
                rns_base,
                glev_context,
            );
        });

        crt_poly.as_mut().copy_from_slice(b_in);
        table.transform_slice(crt_poly.as_mut());
        c_out.neg_assign(self.crt_poly_length, self.poly_length, rns_base.moduli());

        let (_, b_out) = c_out.a_b_mut_slices(self.output_crt_glwe_mid);
        DcrtPolynomial(ArrayBase(b_out)).add_assign(
            &DcrtPolynomial(ArrayBase(b_in)),
            self.poly_length,
            rns_base.moduli(),
        );
    }
}

pub struct CrtGlweKeySwitchingContext<T: UnsignedInteger> {
    big_uint_poly: BigUintPolynomial<Vec<T>>,
    crt_poly: CrtPolynomial<Vec<T>>,
    glev_context: DcrtGlevContext<T>,
}

impl<T: UnsignedInteger> CrtGlweKeySwitchingContext<T> {
    pub fn new(poly_length: usize, crt_poly_len: usize, big_uint_poly_len: usize) -> Self {
        let big_uint_poly = BigUintPolynomial::zero(big_uint_poly_len);
        let crt_poly = CrtPolynomial::zero(crt_poly_len);
        let glev_context = DcrtGlevContext::new(poly_length, crt_poly_len, big_uint_poly_len);
        Self {
            big_uint_poly,
            crt_poly,
            glev_context,
        }
    }

    pub fn as_mut(
        &mut self,
    ) -> (
        &mut BigUintPolynomial<Vec<T>, T>,
        &mut CrtPolynomial<Vec<T>, T>,
        &mut DcrtGlevContext<T>,
    ) {
        (
            &mut self.big_uint_poly,
            &mut self.crt_poly,
            &mut self.glev_context,
        )
    }
}
