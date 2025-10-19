use std::sync::Arc;

use primus_integer::{UnsignedInteger, izip};
use primus_lattice::glev::DcrtGlev;
use primus_modulus::PowOf2Modulus;
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, BigUintPolynomial, Data, DataMut, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;
use primus_reduce::ops::ReduceMul;

use crate::{CrtGlevParameters, CrtGlweCiphertext, CrtGlweSecretKey, DcrtGlweSecretKey};

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
pub struct CrtGlweAutoKey<T, M, Table>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
    Table: DcrtTable<ValueT = T> + Dcrt,
{
    degree: usize,
    poly_length: usize,
    dimension: usize,
    crt_poly_length: usize,
    decompose_length: usize,
    dcrt_glwe_len: usize,
    dcrt_glev_len: usize,
    auto_helper: AutoHelper,
    key: Vec<T>,
    moduli: Vec<M>,
    table: Arc<Table>,
}

impl<T, M, Table> CrtGlweAutoKey<T, M, Table>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
    Table: DcrtTable<ValueT = T> + Dcrt,
{
    pub fn new<R>(
        params: &CrtGlevParameters<T, M>,
        degree: usize,
        sk: &CrtGlweSecretKey<T>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        table: Arc<Table>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
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
                izip!(
                    si.chunks_exact(poly_length),
                    auto_si.chunks_exact_mut(poly_length),
                    moduli
                )
                .for_each(|(a, b, &modulus)| {
                    poly_auto_inplace(a, &auto_helper, b, modulus);
                });

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

                    let mut b_poly = CrtPolynomial(ArrayBase(b));

                    b_poly.add_mul_scalar_residues_assign(
                        &CrtPolynomial(ArrayBase(auto_si.as_ref())),
                        scalar_residues,
                        poly_length,
                        moduli,
                    );

                    let mut b_poly = table.transform_inplace(b_poly);

                    a.chunks_exact_mut(crt_poly_length)
                        .zip(dcrt_sk.iter_dcrt_poly())
                        .for_each(|(ai, si)| {
                            b_poly.add_mul_assign(
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
            decompose_length,
            crt_poly_length,
            dcrt_glwe_len,
            dcrt_glev_len,
            auto_helper,
            key: result,
            moduli: moduli.to_vec(),
            table: Arc::clone(&table),
        }
    }

    pub fn degree(&self) -> usize {
        self.degree
    }

    pub fn table(&self) -> &Table {
        &self.table
    }

    pub fn automorphism_inplace<A, B>(
        &self,
        ciphertext: &CrtGlweCiphertext<A>,
        result: &mut CrtGlweCiphertext<B>,
        params: &CrtGlevParameters<T, M>,
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let crt_poly_length = self.crt_poly_length;
        let dcrt_glwe_len = self.dcrt_glwe_len;
        let dcrt_glev_len = self.dcrt_glev_len;
        let moduli = self.moduli();
        let auto_helper = &self.auto_helper;

        let mut big_uint_poly: BigUintPolynomial<Vec<T>> =
            BigUintPolynomial::zero(poly_length, params.big_uint_value_len());
        let mut temp: CrtGlweCiphertext<Vec<T>> = CrtGlweCiphertext::zero(dcrt_glwe_len);

        izip!(
            temp.iter_crt_poly_mut(crt_poly_length),
            ciphertext.iter_crt_poly(crt_poly_length),
            self.key.chunks_exact(dcrt_glev_len)
        )
        .for_each(|(a, b, key_i)| {
            izip!(
                a.chunks_exact_mut(poly_length),
                b.chunks_exact(poly_length),
                moduli
            )
            .for_each(|(x, y, &modulus)| {
                poly_auto_inplace(y, auto_helper, x, modulus);
            });

            todo!()

            // DcrtGlev::new(ArrayBase(key_i)).mul_polynomial_inplace(
            //     big_uint_polynomial,
            //     result,
            //     dimension,
            //     basis,
            //     table,
            //     rns_base,
            // );
        });

        todo!()
    }

    pub fn moduli_count(&self) -> usize {
        self.moduli.len()
    }

    pub fn moduli(&self) -> &[M] {
        &self.moduli
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
