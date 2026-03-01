use primus_integer::UnsignedInteger;
use primus_lattice::context::DcrtGlevContext;
use primus_lattice::glev::DcrtGlevIterMut;
use primus_ntt::DcrtTable;
use primus_poly::CrtPolynomial;
use primus_reduce::FieldContext;

use crate::{CrtGlevParameters, CrtGlweSecretKey, DcrtGlweSecretKey};

mod coeff;
mod ntt;

pub struct CrtGlweAutoContext<T: UnsignedInteger> {
    crt_poly: CrtPolynomial<Vec<T>>,
    auto_crt_poly: CrtPolynomial<Vec<T>>,
    glev_context: DcrtGlevContext<T>,
}

impl<T: UnsignedInteger> CrtGlweAutoContext<T> {
    pub fn new(poly_length: usize, crt_poly_len: usize, big_uint_poly_len: usize) -> Self {
        let crt_poly = CrtPolynomial::zero(crt_poly_len);
        let auto_crt_poly = CrtPolynomial::zero(crt_poly_len);
        let glev_context = DcrtGlevContext::new(poly_length, crt_poly_len, big_uint_poly_len);
        Self {
            crt_poly,
            auto_crt_poly,
            glev_context,
        }
    }

    pub fn as_mut(
        &mut self,
    ) -> (
        &mut CrtPolynomial<Vec<T>>,
        &mut CrtPolynomial<Vec<T>>,
        &mut DcrtGlevContext<T>,
    ) {
        (
            &mut self.crt_poly,
            &mut self.auto_crt_poly,
            &mut self.glev_context,
        )
    }
}

/// Generate automorphism key data: for each secret-key polynomial s_i,
/// encrypt σ_k(s_i) under a GLEV ciphertext.
///
/// Shared by both [`CrtGlweAutoKey`] (coefficient-domain) and
/// [`DcrtGlweAutoKey`] (NTT-domain) since key generation is identical.
fn generate_auto_key_data<T, M, Table, R>(
    params: &CrtGlevParameters<T, M>,
    coeff_auto_helper: &coeff::CoeffAutoHelper,
    sk: &CrtGlweSecretKey<T>,
    dcrt_sk: &DcrtGlweSecretKey<T>,
    table: &Table,
    rng: &mut R,
) -> Vec<T>
where
    T: UnsignedInteger,
    Table: DcrtTable<ValueT = T>,
    R: rand::Rng + rand::CryptoRng,
    M: FieldContext<T>,
{
    let poly_length = params.poly_length();
    let rns_poly_len = params.rns_poly_len();
    let dcrt_glev_len = params.rns_glev_len();
    let moduli = params.cipher_moduli();

    let mut key = vec![T::ZERO; params.dimension() * dcrt_glev_len];
    let mut auto_si: CrtPolynomial<Vec<T>> = CrtPolynomial::zero(rns_poly_len);

    let key_iter = DcrtGlevIterMut::new(key.as_mut_slice(), dcrt_glev_len);

    sk.iter_crt_poly()
        .zip(key_iter)
        .for_each(|(si, mut dcrt_glev)| {
            coeff::crt_poly_auto_inplace(
                si.0,
                &mut auto_si.0,
                coeff_auto_helper,
                poly_length,
                moduli,
            );

            dcrt_sk.encrypt_dcrt_glev_inplace(&auto_si, &mut dcrt_glev, params, table, rng);
        });

    key
}

pub use coeff::{CoeffAutoHelper, CrtGlweAutoKey, crt_poly_auto_inplace};
pub use ntt::{DcrtGlweAutoKey, NttAutoHelper, dcrt_poly_ntt_auto_inplace};
