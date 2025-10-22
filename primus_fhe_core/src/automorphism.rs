use std::sync::Arc;

use primus_integer::{UnsignedInteger, izip};
use primus_lattice::context::DcrtGlevContext;
use primus_lattice::glev::DcrtGlev;
use primus_modulus::PowOf2Modulus;
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Data, DataMut, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_reduce::ops::ReduceMul;
use primus_rns::RNSBase;

use crate::{
    CrtGlevParameters, CrtGlweCiphertext, CrtGlweSecretKey, DcrtGlweCiphertext, DcrtGlweSecretKey,
};

pub struct CrtGlweAutoContext<T: UnsignedInteger> {
    big_uint_poly: BigUintPolynomial<Vec<T>>,
    auto_crt_poly: CrtPolynomial<Vec<T>>,
    glev_context: DcrtGlevContext<T>,
}

impl<T: UnsignedInteger> CrtGlweAutoContext<T> {
    pub fn new(poly_length: usize, crt_poly_len: usize, big_uint_poly_len: usize) -> Self {
        let big_uint_poly = BigUintPolynomial::zero(big_uint_poly_len);
        let auto_crt_poly = CrtPolynomial::zero(crt_poly_len);
        let glev_context = DcrtGlevContext::new(poly_length, crt_poly_len, big_uint_poly_len);
        Self {
            big_uint_poly,
            auto_crt_poly,
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
    DimensionPlusOne,
    One,
}

/// Automorphism key
#[derive(Clone)]
pub struct CrtGlweAutoKey<T, Table>
where
    T: UnsignedInteger,
    Table: DcrtTable<ValueT = T> + Dcrt,
{
    degree: usize,
    poly_length: usize,
    dimension: usize,
    crt_poly_length: usize,
    dcrt_glwe_mid: usize,
    dcrt_glwe_len: usize,
    dcrt_glev_len: usize,
    auto_helper: AutoHelper,
    key: Vec<T>,
    table: Arc<Table>,
}

impl<T, Table> CrtGlweAutoKey<T, Table>
where
    T: UnsignedInteger,
    Table: DcrtTable<ValueT = T> + Dcrt,
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
        let poly_length = sk.poly_length();
        let dimension = sk.dimension();
        let moduli_count = sk.moduli_count();
        let crt_poly_length = sk.crt_poly_length();
        let dcrt_glwe_len = dcrt_sk.crt_glwe_len();
        let basis = params.basis();
        let decompose_length = basis.decompose_length();
        let dcrt_glev_len = decompose_length * dcrt_glwe_len;
        let moduli = params.cipher_moduli();
        let moduli_value = params.cipher_moduli_value();
        let gaussian = params.noise_distribution();
        let uniform_distrs = params.cipher_moduli_uniform_distr();

        let dcrt_glwe_mid = crt_poly_length * dimension;

        debug_assert_eq!(moduli_count, params.cipher_moduli().len());

        let auto_helper = if degree == 1 {
            AutoHelper::One
        } else if degree == poly_length + 1 {
            AutoHelper::DimensionPlusOne
        } else {
            AutoHelper::Permutation(generate_permutate_ops(degree, poly_length))
        };

        let mut result = vec![T::ZERO; dimension * dcrt_glev_len];
        let mut auto_si = vec![T::ZERO; crt_poly_length];

        result
            .chunks_exact_mut(dcrt_glev_len)
            .zip(sk.iter_crt_poly())
            .for_each(|(dcrt_glev, si)| {
                crt_poly_auto_inplace(si, &mut auto_si, &auto_helper, poly_length, moduli);

                izip!(
                    dcrt_glev.chunks_exact_mut(dcrt_glwe_len),
                    basis.scalars_residues_iter()
                )
                .for_each(|(dcrt_glwe, scalar_residues)| {
                    let (a, b) = unsafe { dcrt_glwe.split_at_mut_unchecked(dcrt_glwe_mid) };
                    primus_distr::sample_crt_gaussian_values_inplace(
                        b,
                        poly_length,
                        moduli_value,
                        gaussian,
                        rng,
                    );

                    let mut b_crt_poly = CrtPolynomial(ArrayBase(b));

                    b_crt_poly.add_mul_scalar_assign(
                        &CrtPolynomial(ArrayBase(auto_si.as_ref())),
                        scalar_residues,
                        poly_length,
                        moduli,
                    );

                    let mut b_dcrt_poly = table.transform_inplace(b_crt_poly);

                    a.chunks_exact_mut(crt_poly_length)
                        .zip(dcrt_sk.iter_dcrt_poly())
                        .for_each(|(ai, si)| {
                            primus_distr::sample_crt_uniform_values_inplace(
                                ai,
                                poly_length,
                                uniform_distrs,
                                rng,
                            );
                            b_dcrt_poly.add_mul_assign(
                                &DcrtPolynomial(ArrayBase(ai)),
                                &DcrtPolynomial(ArrayBase(si)),
                                poly_length,
                                moduli,
                            );
                        });
                });
            });

        Self {
            degree,
            poly_length,
            dimension,
            crt_poly_length,
            dcrt_glwe_len,
            dcrt_glwe_mid,
            dcrt_glev_len,
            auto_helper,
            key: result,
            table: Arc::clone(&table),
        }
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    pub fn iter_dcrt_glev(&self) -> std::slice::ChunksExact<'_, T> {
        self.key.chunks_exact(self.dcrt_glev_len)
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
        let poly_length = self.poly_length;
        let crt_poly_length = self.crt_poly_length;
        let dcrt_glwe_mid = self.dcrt_glwe_mid;
        let moduli = params.cipher_moduli();
        let auto_helper = &self.auto_helper;

        let (big_uint_poly, auto_crt_poly, glev_context) = context.as_mut();

        result.set_zero();
        let mut temp = DcrtGlweCiphertext::new(ArrayBase(result.as_mut()));

        let (a_in, b_in) = ciphertext.a_b_slices(dcrt_glwe_mid);

        izip!(a_in.chunks_exact(crt_poly_length), self.iter_dcrt_glev()).for_each(
            |(in_crt_poly, auto_key_i)| {
                crt_poly_auto_inplace(
                    in_crt_poly,
                    auto_crt_poly.as_mut(),
                    auto_helper,
                    poly_length,
                    moduli,
                );

                rns_base.compose_polynomial_inplace(&auto_crt_poly, big_uint_poly, poly_length);

                let auto_key = DcrtGlev::new(ArrayBase(auto_key_i));

                temp.add_dcrt_glev_mul_big_uint_poly_assign(
                    &auto_key,
                    &big_uint_poly,
                    params.basis(),
                    self.table(),
                    rns_base,
                    glev_context,
                );
            },
        );

        crt_poly_auto_inplace(
            b_in,
            auto_crt_poly.as_mut(),
            auto_helper,
            poly_length,
            moduli,
        );

        let (a_out, b_out) = result.a_b_mut_slices(dcrt_glwe_mid);

        a_out
            .chunks_exact_mut(crt_poly_length)
            .for_each(|crt_poly| {
                let mut temp = CrtPolynomial(ArrayBase(crt_poly));
                temp.neg_assign(poly_length, moduli);
                self.table.transform_inplace(temp);
            });

        let mut temp = CrtPolynomial(ArrayBase(b_out));
        temp.sub_assign(auto_crt_poly, poly_length, moduli);
        self.table.transform_inplace(temp);
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

fn crt_poly_auto_inplace<T, M>(
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
    .for_each(|(in_poly, auto_poly, &modulus)| {
        poly_auto_inplace(in_poly, auto_helper, auto_poly, modulus);
    });
}

#[inline]
fn poly_auto_inplace<T, M>(poly: &[T], auto_helper: &AutoHelper, result: &mut [T], modulus: M)
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    match auto_helper {
        AutoHelper::Permutation(from_ops) => {
            poly_auto_inplace_for_permutation(poly, result, from_ops, modulus);
        }
        AutoHelper::DimensionPlusOne => {
            poly_auto_inplace_for_dimension_plus_one(poly, result, modulus);
        }
        AutoHelper::One => poly_auto_inplace_for_one(poly, result),
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
