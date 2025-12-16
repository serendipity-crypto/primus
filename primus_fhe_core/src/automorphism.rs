use std::sync::Arc;

use primus_integer::{Data, DataMut, RawData, UnsignedInteger, izip};
use primus_lattice::glev::DcrtGlevIter;
use primus_lattice::{context::DcrtGlevContext, glev::DcrtGlevIterMut};
use primus_modulus::PowOf2Modulus;
use primus_ntt::DcrtTable;
use primus_poly::{CrtPolynomial, DcrtPolynomial};
use primus_reduce::FieldContext;
use primus_reduce::ops::ReduceMul;
use primus_rns::RNSBase;

use crate::{
    CrtGlevParameters, CrtGlweCiphertext, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweSecretKey,
};

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

/// This defines the operation when perform automorphism on each coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Add,
    Sub,
}

/// This defines the operation and the source index
/// when perform automorphism on each coefficient.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FromOp {
    from: usize,
    op: Op,
}

#[derive(Debug, Clone)]
pub enum AutoHelper {
    Permutation(Vec<FromOp>),
    PolyLengthPlusOne,
    One,
}

/// Automorphism key
#[derive(Clone)]
pub struct CrtGlweAutoKey<T, Table>
where
    T: UnsignedInteger,
    Table: DcrtTable<ValueT = T>,
{
    key: Vec<T>,
    degree: usize,
    rns_glev_len: usize,
    auto_helper: AutoHelper,
    table: Arc<Table>,
}

impl<T, Table> CrtGlweAutoKey<T, Table>
where
    T: UnsignedInteger,
    Table: DcrtTable<ValueT = T>,
{
    pub fn new<M, R>(
        params: &CrtGlevParameters<T, M>,
        degree: usize,
        sk: &CrtGlweSecretKey<T>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        table: Arc<Table>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        M: FieldContext<T>,
    {
        let poly_length = params.poly_length();
        let rns_poly_len = params.rns_poly_len();
        let dcrt_glev_len = params.rns_glev_len();
        let moduli = params.cipher_moduli();

        let auto_helper = if degree == 1 {
            AutoHelper::One
        } else if degree == poly_length + 1 {
            AutoHelper::PolyLengthPlusOne
        } else {
            AutoHelper::Permutation(generate_permutate_ops(degree, poly_length))
        };

        let mut key = vec![T::ZERO; params.dimension() * dcrt_glev_len];
        let mut auto_si: CrtPolynomial<Vec<T>> = CrtPolynomial::zero(rns_poly_len);

        let key_iter = DcrtGlevIterMut::new(key.as_mut_slice(), dcrt_glev_len);

        sk.iter_crt_poly()
            .zip(key_iter)
            .for_each(|(si, mut dcrt_glev)| {
                crt_poly_auto_inplace(si.0, &mut auto_si.0, &auto_helper, poly_length, moduli);

                dcrt_sk.encrypt_dcrt_glev_inplace(
                    &auto_si,
                    &mut dcrt_glev,
                    params,
                    table.as_ref(),
                    rng,
                );
            });

        Self {
            key,
            degree,
            rns_glev_len: dcrt_glev_len,
            auto_helper,
            table: Arc::clone(&table),
        }
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn auto_helper(&self) -> &AutoHelper {
        &self.auto_helper
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    pub fn iter_dcrt_glev(&self) -> DcrtGlevIter<'_, T> {
        DcrtGlevIter::new(self.key.as_slice(), self.rns_glev_len)
    }

    pub fn automorphism_inplace<M, A, B>(
        &self,
        ciphertext: &CrtGlweCiphertext<A>,
        result: &mut CrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweAutoContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let rns_glwe_mid = params.rns_glwe_mid();
        let moduli = params.cipher_moduli();

        let auto_helper = &self.auto_helper;

        debug_assert_eq!(ciphertext.as_ref().len(), params.rns_glwe_len());

        let (_, auto_crt_poly, glev_context) = context.as_mut();

        result.set_zero();
        let mut temp = DcrtGlweCiphertext::new(result.as_mut());

        let (a_in, b_in) = ciphertext.a_b(rns_glwe_mid);

        self.iter_dcrt_glev()
            .zip(a_in)
            .for_each(|(auto_key_i, in_crt_poly)| {
                crt_poly_auto_inplace(
                    in_crt_poly.0,
                    auto_crt_poly.as_mut(),
                    auto_helper,
                    poly_length,
                    moduli,
                );

                temp.add_dcrt_glev_mul_crt_poly_assign(
                    &auto_key_i,
                    auto_crt_poly,
                    params.basis(),
                    self.table(),
                    rns_base,
                    glev_context,
                );
            });

        crt_poly_auto_inplace(
            b_in.0,
            auto_crt_poly.as_mut(),
            auto_helper,
            poly_length,
            moduli,
        );

        let _ = temp.into_coeff_form(self.table());

        let (a_out, mut b_out) = result.a_b_mut(rns_glwe_mid);

        a_out.for_each(|mut ai| ai.neg_assign(poly_length, moduli));

        auto_crt_poly.sub_to_right(&mut b_out, poly_length, moduli);
    }

    pub fn automorphism_to_dcrt_glwe_inplace<M, A, B>(
        &self,
        ciphertext: &DcrtGlweCiphertext<A>,
        result: &mut DcrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
        rns_base: &RNSBase<T, M>,
        context: &mut CrtGlweAutoContext<T>,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = params.poly_length();
        let rns_glwe_mid = params.rns_glwe_mid();
        let moduli = params.cipher_moduli();

        let auto_helper = &self.auto_helper;

        debug_assert_eq!(ciphertext.as_ref().len(), params.rns_glwe_len());

        let (crt_poly, auto_crt_poly, glev_context) = context.as_mut();

        result.set_zero();

        let (a_in, b_in) = ciphertext.a_b(rns_glwe_mid);

        self.iter_dcrt_glev()
            .zip(a_in)
            .for_each(|(auto_key_i, in_dcrt_poly)| {
                crt_poly.as_mut().copy_from_slice(in_dcrt_poly.0);
                self.table.inverse_transform_slice(crt_poly.as_mut());

                crt_poly_auto_inplace(
                    crt_poly.as_ref(),
                    auto_crt_poly.as_mut(),
                    auto_helper,
                    poly_length,
                    moduli,
                );

                result.add_dcrt_glev_mul_crt_poly_assign(
                    &auto_key_i,
                    auto_crt_poly,
                    params.basis(),
                    self.table(),
                    rns_base,
                    glev_context,
                );
            });

        crt_poly.as_mut().copy_from_slice(b_in.0);
        self.table.inverse_transform_slice(crt_poly.as_mut());

        crt_poly_auto_inplace(
            crt_poly.as_ref(),
            auto_crt_poly.as_mut(),
            auto_helper,
            poly_length,
            moduli,
        );

        self.table.transform_slice(auto_crt_poly.as_mut());

        let (a_out, mut b_out) = result.a_b_mut(rns_glwe_mid);

        a_out.for_each(|mut ai| ai.neg_assign(poly_length, moduli));

        DcrtPolynomial(auto_crt_poly.as_ref()).sub_to_right(&mut b_out, poly_length, moduli);
    }
}

#[inline]
fn generate_permutate_ops(degree: usize, poly_length: usize) -> Vec<FromOp> {
    let twice_poly_length = poly_length << 1;
    let modulus = <PowOf2Modulus<usize>>::new(twice_poly_length);

    let mut result = vec![
        FromOp {
            from: 0,
            op: Op::Add
        };
        poly_length
    ];

    for i in 0..poly_length {
        let to = modulus.reduce_mul(i, degree);
        if to < poly_length {
            result[to] = FromOp {
                from: i,
                op: Op::Add,
            };
        } else {
            result[to - poly_length] = FromOp {
                from: i,
                op: Op::Sub,
            };
        }
    }
    result
}

pub fn crt_poly_auto_inplace<T, M>(
    crt_poly: &[T],
    auto_crt_poly: &mut [T],
    auto_helper: &AutoHelper,
    poly_length: usize,
    moduli: &[M],
) where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    izip!(
        crt_poly.chunks_exact(poly_length),
        auto_crt_poly.chunks_exact_mut(poly_length),
        moduli
    )
    .for_each(|(poly, auto_poly, &modulus)| {
        poly_auto_inplace(poly, auto_poly, auto_helper, modulus);
    });
}

#[inline]
fn poly_auto_inplace<T, M>(poly: &[T], auto_poly: &mut [T], auto_helper: &AutoHelper, modulus: M)
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    match auto_helper {
        AutoHelper::Permutation(from_ops) => {
            poly_auto_inplace_for_permutation(poly, auto_poly, from_ops, modulus);
        }
        AutoHelper::PolyLengthPlusOne => {
            poly_auto_inplace_for_dimension_plus_one(poly, auto_poly, modulus);
        }
        AutoHelper::One => poly_auto_inplace_for_one(poly, auto_poly),
    }
}

#[inline]
fn poly_auto_inplace_for_permutation<T, M>(
    poly: &[T],
    result: &mut [T],
    from_ops: &[FromOp],
    modulus: M,
) where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    for (d, from_op) in result.iter_mut().zip(from_ops.iter()) {
        let c = unsafe { poly.get_unchecked(from_op.from) };
        match from_op.op {
            Op::Add => {
                *d = *c;
            }
            Op::Sub => {
                *d = modulus.reduce_neg(*c);
            }
        }
    }
}

#[inline]
fn poly_auto_inplace_for_dimension_plus_one<T, M>(poly: &[T], result: &mut [T], modulus: M)
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    for (pi, di) in unsafe {
        poly.as_chunks_unchecked::<2>()
            .iter()
            .zip(result.as_chunks_unchecked_mut::<2>())
    } {
        di[0] = pi[0];
        di[1] = modulus.reduce_neg(pi[1]);
    }
}

#[inline]
fn poly_auto_inplace_for_one<T>(poly: &[T], result: &mut [T])
where
    T: UnsignedInteger,
{
    result.copy_from_slice(poly);
}
