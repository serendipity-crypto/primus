use std::sync::Arc;

use primus_decompose::big_integer::BigUintApproxSignedBasis;
use primus_distr::SignedDiscreteGaussian;
use primus_integer::{UnsignedInteger, izip};
use primus_modulus::PowOf2Modulus;
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, NttPolynomial, RawData};
use primus_reduce::FieldContext;
use primus_reduce::ops::ReduceMul;
use rand::distr::Distribution;

use crate::{CrtGlweCiphertext, CrtGlweSecretKey, DcrtGlweSecretKey};

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
    decompose_length: usize,
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
        sk: &CrtGlweSecretKey<T>,
        dcrt_sk: &DcrtGlweSecretKey<T>,
        degree: usize,
        gaussian: &SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
        basis: &BigUintApproxSignedBasis<T>,
        moduli: &[M],
        table: Arc<Table>,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
    {
        let poly_length = sk.poly_length();
        let dimension = sk.dimension();
        let moduli_count = sk.moduli_count();
        let decompose_length = basis.decompose_length();

        debug_assert_eq!(moduli_count, moduli.len());

        let a_b_mid = dimension * poly_length;
        let glwe_len = a_b_mid + poly_length;
        let glev_len = decompose_length * glwe_len;
        let single_modulus_len = dimension * glev_len;

        let e_glev_len = decompose_length * poly_length;
        let e_single_modulus_len = dimension * e_glev_len;

        let auto_helper = if degree == 1 {
            AutoHelper::One
        } else if degree == poly_length + 1 {
            AutoHelper::DimensionPlusOne
        } else {
            AutoHelper::Permutation(generate_permutate_ops(degree, poly_length))
        };

        let mut result = vec![T::ZERO; moduli_count * single_modulus_len];
        let mut e_all = vec![T::ZERO; moduli_count * e_single_modulus_len];
        let modulus_values: Vec<T> = moduli.iter().map(|m| m.value_unchecked()).collect();

        primus_distr::sample_crt_gaussian_values_inplace(
            &mut e_all,
            e_single_modulus_len,
            &modulus_values,
            gaussian,
            rng,
        );

        izip!(
            result.chunks_exact_mut(single_modulus_len),
            sk.iter_each_modulus(),
            dcrt_sk.iter_each_modulus(),
            e_all.chunks_exact_mut(e_single_modulus_len),
            basis.scalars_residue().chunks_exact(decompose_length),
            table.iter(),
            moduli,
        )
        .for_each(|(auto_key, key, ntt_key, es, sclars, ntt_table, modulus)| {
            let uniform_distr = modulus.uniform_distribution();
            izip!(
                auto_key.chunks_exact_mut(glev_len),
                key.chunks_exact(poly_length),
                es.chunks_exact_mut(e_glev_len)
            )
            .for_each(|(glev, key_part, e_glev)| {
                izip!(
                    glev.chunks_exact_mut(glwe_len),
                    e_glev.chunks_exact_mut(poly_length),
                    sclars
                )
                .for_each(|(glwe, e_glwe, scalar)| {
                    let (a, b) = unsafe { glwe.split_at_mut_unchecked(a_b_mid) };

                    b.copy_from_slice(e_glwe);

                    poly_auto_inplace(key_part, &auto_helper, e_glwe, *modulus);
                    ArrayBase(&mut *b).add_mul_scalar_assign(&ArrayBase(e_glwe), *scalar, *modulus);
                    ntt_table.transform_slice(b);

                    a.iter_mut()
                        .zip(uniform_distr.sample_iter(&mut *rng))
                        .for_each(|(i, o)| *i = o);

                    let mut b_poly = NttPolynomial(ArrayBase(b));

                    a.chunks_exact_mut(poly_length)
                        .zip(ntt_key.chunks_exact(poly_length))
                        .for_each(|(ai, s)| {
                            b_poly.add_mul_assign(
                                &NttPolynomial(ArrayBase(ai)),
                                &NttPolynomial(ArrayBase(s)),
                                *modulus,
                            );
                        });
                });
            });
        });

        Self {
            degree,
            poly_length,
            dimension,
            decompose_length,
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
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = self.poly_length;
        let single_modulus_len = self.dimension * poly_length;

        izip!(
            result.iter_each_modulus_mut(single_modulus_len),
            ciphertext.iter_each_modulus(single_modulus_len),
            self.moduli(),
        )
        .for_each(|(glwe_out, glwe_in, modulus)| {
            glwe_out
                .chunks_exact_mut(poly_length)
                .zip(glwe_in.chunks_exact(poly_length))
                .for_each(|(poly_out, poly_in)| {
                    poly_auto_inplace(poly_in, &self.auto_helper, poly_out, *modulus);
                });
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
