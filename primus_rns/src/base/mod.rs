use std::slice::Iter;

use itertools::Itertools;
use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{
    BigIntegerOps, UnsignedInteger, izip, multiply_many_values, multiply_many_values_except_inplace,
};
use primus_modulo::ops::*;
use primus_poly::{BigUintPolynomial, Data, DataMut, RawData, crt::CrtPolynomial};
use primus_reduce::FieldContext;

use crate::RNSError;

/// A residue number system or residue numeral system (RNS) is a numeral system representing integers
/// by their values modulo several pairwise coprime integers called the moduli.
/// This representation is allowed by the Chinese remainder theorem,
/// which asserts that, if M is the product of the moduli, there is,
/// in an interval of length M, exactly one integer having any given set of modular values.
/// Using a residue numeral system for arithmetic operations is also called multi-modular arithmetic.
#[derive(Clone)]
pub struct RNSBase<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    pub moduli: Vec<M>,
    pub moduli_product: Vec<T>,
    pub punctured_product: Vec<T>,
    pub inv_punctured_product_mod_modulus: Vec<ShoupFactor<T>>,
}

impl<T, M> RNSBase<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Creates a new [`RNSBase<T, M>`].
    ///
    /// # Panics
    ///
    /// Panics if any inverse modulo operation panics.
    ///
    /// # Errors
    ///
    /// This function will return an error if moduli are not co-prime with each others.
    #[inline]
    pub fn new(moduli: &[M]) -> Result<Self, RNSError> {
        let modulus_values = moduli
            .iter()
            .map(|m| m.value_unchecked())
            .collect::<Vec<_>>();

        if modulus_values
            .iter()
            .tuple_combinations()
            .any(|(&a, &b)| a.not_coprime(b))
        {
            return Err(RNSError::CoPrimeError);
        }

        let moduli_product = multiply_many_values(&modulus_values);

        let big_uint_len = moduli_product.len();
        let mut punctured_product = vec![T::ZERO; big_uint_len * moduli.len()];
        punctured_product
            .chunks_exact_mut(big_uint_len)
            .enumerate()
            .for_each(|(i, chunk)| {
                multiply_many_values_except_inplace(&modulus_values, i, chunk);
            });

        let inv_punctured_product_mod_modulus = punctured_product
            .chunks_exact(big_uint_len)
            .zip(moduli)
            .map(|(p, &modulus)| {
                let inv = p.modulo(modulus).try_inv_modulo(modulus).unwrap();
                ShoupFactor::new(inv, modulus.value_unchecked())
            })
            .collect::<Vec<ShoupFactor<T>>>();

        Ok(Self {
            moduli: moduli.to_vec(),
            moduli_product,
            punctured_product,
            inv_punctured_product_mod_modulus,
        })
    }

    /// Returns a reference to the moduli of this [`RNSBase<T, M>`].
    #[inline]
    pub fn moduli(&self) -> &[M] {
        &self.moduli
    }

    /// Returns a reference to the moduli product of this [`RNSBase<T, M>`].
    #[inline]
    pub fn moduli_product(&self) -> &[T] {
        &self.moduli_product
    }

    /// Returns a reference to the punctured product of this [`RNSBase<T, M>`].
    #[inline]
    pub fn punctured_product(&self) -> &[T] {
        &self.punctured_product
    }

    /// Returns an iterator over the punctured product of this [`RNSBase<T, M>`].
    #[inline]
    pub fn punctured_product_iter(&self) -> std::slice::ChunksExact<'_, T> {
        self.punctured_product
            .chunks_exact(self.moduli_product.len())
    }

    /// Returns a reference to the inverse punctured product mod modulus of this [`RNSBase<T, M>`].
    #[inline]
    pub fn inv_punctured_product_mod_modulus(&self) -> &[ShoupFactor<T>] {
        &self.inv_punctured_product_mod_modulus
    }

    /// Decomposes a value into its RNS representation.
    #[inline]
    pub fn decompose(&self, value: &[T]) -> Vec<T> {
        self.moduli
            .iter()
            .map(|&modulus| value.modulo(modulus))
            .collect()
    }

    /// Decomposes a value into its RNS representation, writing the result into the provided slice.
    #[inline]
    pub fn decompose_inplace(&self, value: &[T], residues: &mut [T]) {
        debug_assert_eq!(self.moduli.len(), residues.len());

        for (residue, &modulus) in residues.iter_mut().zip(self.moduli.iter()) {
            *residue = value.modulo(modulus);
        }
    }

    pub fn decompose_multiple_values_inplace(
        &self,
        big_uint_values: &[T],
        multi_residues: &mut [T],
        value_count: usize,
    ) {
        let value_len = self.moduli_product.len();
        for (residues, &modulus) in multi_residues
            .chunks_exact_mut(value_count)
            .zip(self.moduli())
        {
            for (residue, value) in residues
                .iter_mut()
                .zip(big_uint_values.chunks_exact(value_len))
            {
                *residue = value.modulo(modulus);
            }
        }
    }

    pub fn decompose_polynomial_inplace<R, W>(
        &self,
        big_uint_poly: &BigUintPolynomial<R>,
        crt_poly: &mut CrtPolynomial<W>,
        poly_length: usize,
    ) where
        R: RawData<Elem = T> + Data,
        W: RawData<Elem = T> + DataMut,
    {
        let value_len = self.moduli_product.len();
        for (poly, &modulus) in crt_poly
            .iter_each_modulus_mut(poly_length)
            .zip(self.moduli())
        {
            for (residue, value) in poly.iter_mut().zip(big_uint_poly.iter(value_len)) {
                *residue = value.modulo(modulus);
            }
        }
    }

    /// Composes a value from its RNS representation.
    pub fn compose(&self, residues: &[T]) -> Vec<T> {
        debug_assert_eq!(self.moduli.len(), residues.len());

        let value_len = self.moduli_product.len();

        let mut value = vec![T::ZERO; value_len];

        izip!(
            residues,
            &self.inv_punctured_product_mod_modulus,
            self.punctured_product.chunks_exact(value_len),
            &self.moduli
        )
        .for_each(
            |(&ri, &inv_mi, mi, &modulus): (&T, &ShoupFactor<T>, &[T], &M)| {
                let product = inv_mi.factor_mul_modulo(ri, modulus.value_unchecked());
                let carry = mi.slice_mul_value_add_inplace(product, &mut value);
                if !carry.is_zero() || value.slice_cmp(&self.moduli_product).is_ge() {
                    let _ = value.slice_sub_assign(&self.moduli_product);
                }
            },
        );

        value
    }

    pub fn compose_inplace(&self, residues: &[T], value: &mut [T]) {
        debug_assert_eq!(self.moduli.len(), residues.len());
        debug_assert_eq!(self.moduli_product.len(), value.len());

        let value_len = self.moduli_product.len();

        izip!(
            residues,
            &self.inv_punctured_product_mod_modulus,
            self.punctured_product.chunks_exact(value_len),
            &self.moduli
        )
        .for_each(
            |(&ri, &inv_mi, mi, &modulus): (&T, &ShoupFactor<T>, &[T], &M)| {
                let product = inv_mi.factor_mul_modulo(ri, modulus.value_unchecked());
                let carry = mi.slice_mul_value_add_inplace(product, value);
                if !carry.is_zero() || value.slice_cmp(&self.moduli_product).is_ge() {
                    let _ = value.slice_sub_assign(&self.moduli_product);
                }
            },
        );
    }

    pub fn compose_multiple_values_inplace(
        &self,
        multi_residues: &[T],
        big_uint_values: &mut [T],
        value_count: usize,
    ) {
        let value_len = self.moduli_product.len();
        let mut residues = vec![T::ZERO; self.moduli.len()];

        let mut iters: Vec<Iter<'_, T>> = multi_residues
            .chunks_exact(value_count)
            .map(|s| s.iter())
            .collect();

        for value in big_uint_values.chunks_exact_mut(value_len) {
            for (iter, residue) in iters.iter_mut().zip(residues.iter_mut()) {
                *residue = *iter.next().unwrap();
            }
            self.compose_inplace(&residues, value);
        }
    }

    pub fn compose_polynomial_inplace<R, W>(
        &self,
        crt_poly: &CrtPolynomial<R>,
        big_uint_poly: &mut BigUintPolynomial<W>,
        poly_length: usize,
    ) where
        R: RawData<Elem = T> + Data,
        W: RawData<Elem = T> + DataMut,
    {
        let value_len = self.moduli_product.len();

        let mut residues = vec![T::ZERO; self.moduli.len()];
        let mut iters: Vec<Iter<'_, T>> = crt_poly
            .iter_each_modulus(poly_length)
            .map(|s| s.iter())
            .collect();
        for value in big_uint_poly.iter_mut(value_len) {
            for (iter, res) in iters.iter_mut().zip(residues.iter_mut()) {
                *res = *iter.next().unwrap();
            }
            self.compose_inplace(&residues, value);
        }
    }
}
