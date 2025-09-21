use std::{
    ops::{Index, IndexMut},
    slice::SliceIndex,
};

use integer::{
    BigIntegerOps, UnsignedInteger, izip, multiply_many_values, multiply_many_values_except_inplace,
};
use itertools::Itertools;
use modulo::ops::*;
use primus_factor::{FactorMul, ShoupFactor};
use reduce::FieldAdapter;

use crate::RNSError;

#[derive(Clone)]
pub struct RNSBase<T: UnsignedInteger, M: FieldAdapter<T>> {
    pub base: Vec<M>,
    pub base_product: Vec<T>,
    pub punctured_product: Vec<T>,
    pub inv_punctured_product_mod_base: Vec<ShoupFactor<T>>,
}

impl<T, M, I> Index<I> for RNSBase<T, M>
where
    T: UnsignedInteger,
    M: FieldAdapter<T>,
    I: SliceIndex<[M]>,
{
    type Output = I::Output;

    #[inline]
    fn index(&self, index: I) -> &Self::Output {
        Index::index(&*self.base, index)
    }
}

impl<T, M, I> IndexMut<I> for RNSBase<T, M>
where
    T: UnsignedInteger,
    M: FieldAdapter<T>,
    I: SliceIndex<[M]>,
{
    #[inline]
    fn index_mut(&mut self, index: I) -> &mut Self::Output {
        IndexMut::index_mut(&mut *self.base, index)
    }
}

impl<T: UnsignedInteger, M: FieldAdapter<T>> RNSBase<T, M> {
    /// Create a new RNS base from the given moduli.
    pub fn new(base: &[M]) -> Result<Self, RNSError> {
        let base_values = base.iter().map(|m| m.value_unchecked()).collect::<Vec<_>>();

        if base_values
            .iter()
            .tuple_combinations()
            .any(|(&a, &b)| a.not_coprime(b))
        {
            return Err(RNSError::CoPrimeError);
        }

        let base_product = multiply_many_values(&base_values);

        let chunk_size = base_product.len();
        let len = chunk_size * base.len();
        let mut punctured_product = vec![T::ZERO; len];
        punctured_product
            .chunks_exact_mut(chunk_size)
            .enumerate()
            .for_each(|(i, chunk)| {
                multiply_many_values_except_inplace(&base_values, i, chunk);
            });

        let inv_punctured_product_mod_base = punctured_product
            .chunks_exact(chunk_size)
            .zip(base.iter())
            .map(|(p, &m)| {
                let inv = p.modulo(m).try_inv_modulo(m).unwrap();
                ShoupFactor::new(inv, m.value_unchecked())
            })
            .collect::<Vec<ShoupFactor<T>>>();

        Ok(Self {
            base: base.to_vec(),
            base_product,
            punctured_product,
            inv_punctured_product_mod_base,
        })
    }

    pub fn base(&self) -> &[M] {
        &self.base
    }

    pub fn base_product(&self) -> &[T] {
        &self.base_product
    }

    pub fn punctured_product(&self) -> &[T] {
        &self.punctured_product
    }

    pub fn punctured_product_iter(&self) -> std::slice::ChunksExact<'_, T> {
        self.punctured_product.chunks_exact(self.base_product.len())
    }

    pub fn inv_punctured_product_mod_base(&self) -> &[ShoupFactor<T>] {
        &self.inv_punctured_product_mod_base
    }

    /// Decomposes a value into its RNS representation.
    pub fn decompose(&self, value: &[T]) -> Vec<T> {
        self.base.iter().map(|&m| value.modulo(m)).collect()
    }

    /// Decomposes a value into its RNS representation, writing the result into the provided slice.
    pub fn decompose_inplace(&self, value: &[T], residues: &mut [T]) {
        debug_assert_eq!(self.base.len(), residues.len());

        for (r, &m) in residues.iter_mut().zip(self.base.iter()) {
            *r = value.modulo(m);
        }
    }

    /// Composes a value from its RNS representation.
    pub fn compose(&self, residues: &[T]) -> Vec<T> {
        debug_assert_eq!(self.base.len(), residues.len());

        let mut value = vec![T::ZERO; self.base_product.len()];

        izip!(
            residues,
            &self.inv_punctured_product_mod_base,
            self.punctured_product.chunks_exact(self.base_product.len()),
            &self.base
        )
        .for_each(
            |(&ri, &inv_mi, mi, &modulus): (&T, &ShoupFactor<T>, &[T], &M)| {
                let product = inv_mi.factor_mul_modulo(ri, modulus.value_unchecked());
                let carry = mi.slice_mul_value_add_inplace(product, &mut value);
                if !carry.is_zero() || value.slice_cmp(&self.base_product).is_ge() {
                    let _ = value.slice_sub_assign(&self.base_product);
                }
            },
        );

        value
    }

    pub fn compose_inplace(&self, residues: &[T], value: &mut [T]) {
        debug_assert_eq!(self.base.len(), residues.len());
        debug_assert_eq!(self.base.len(), value.len());

        izip!(
            residues,
            &self.inv_punctured_product_mod_base,
            self.punctured_product.chunks_exact(self.base_product.len()),
            &self.base
        )
        .for_each(
            |(&ri, &inv_mi, mi, &modulus): (&T, &ShoupFactor<T>, &[T], &M)| {
                let product = inv_mi.factor_mul_modulo(ri, modulus.value_unchecked());
                let carry = mi.slice_mul_value_add_inplace(product, value);
                if !carry.is_zero() || value.slice_cmp(&self.base_product).is_ge() {
                    let _ = value.slice_sub_assign(&self.base_product);
                }
            },
        );
    }
}
