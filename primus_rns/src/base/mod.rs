use crate::slice::BigIntOps;
use integer::UnsignedInteger;
use integer::izip;
use itertools::Itertools;
use modulo::ops::*;
use primus_factor::FactorMul;
use primus_factor::ShoupFactor;
use reduce::FieldAdapter;

use crate::{
    RNSError,
    slice::{multiply_many_values, multiply_many_values_except},
};

pub struct RNSBase<T: UnsignedInteger, M: FieldAdapter<T>> {
    pub moduli: Vec<M>,
    pub moduli_product: Vec<T>,
    pub punctured_moduli: Vec<Vec<T>>,
    pub inv_punctured_moduli: Vec<ShoupFactor<T>>,
}

impl<T: UnsignedInteger, M: FieldAdapter<T>> RNSBase<T, M> {
    pub fn new(moduli: &[M]) -> Result<Self, RNSError> {
        let moduli_values = moduli
            .iter()
            .map(|m| m.value_unchecked())
            .collect::<Vec<_>>();

        if moduli_values
            .iter()
            .tuple_combinations()
            .any(|(&a, &b)| a.not_coprime(b))
        {
            return Err(RNSError::CoPrimeError);
        }

        let moduli_product = multiply_many_values(&moduli_values);
        let punctured_moduli = (0..moduli.len())
            .map(|i| multiply_many_values_except(&moduli_values, i))
            .collect::<Vec<_>>();
        let inv_punctured_moduli = punctured_moduli
            .iter()
            .zip(moduli.iter())
            .map(|(p, &m)| {
                let inv = p.as_slice().modulo(m).try_inv_modulo(m).unwrap();
                ShoupFactor::new(inv, m.value_unchecked())
            })
            .collect::<Vec<ShoupFactor<T>>>();

        Ok(Self {
            moduli: moduli.to_vec(),
            moduli_product,
            punctured_moduli,
            inv_punctured_moduli,
        })
    }

    pub fn decompose(&self, value: &[T]) -> Vec<T> {
        let mut result = vec![T::ZERO; self.moduli.len()];

        if self.moduli.len() > 1 {
            for (r, &b) in result.iter_mut().zip(self.moduli.iter()) {
                *r = value.modulo(b);
            }
        } else {
            result[0] = value[0];
        }

        result
    }

    pub fn decompose_inplace(&self, value: &[T], result: &mut [T]) {
        if self.moduli.len() > 1 {
            for (r, &b) in result.iter_mut().zip(self.moduli.iter()) {
                *r = value.modulo(b);
            }
        } else {
            result[0] = value[0];
        }
    }

    pub fn compose(&self, residues: &[T]) -> Vec<T> {
        assert_eq!(residues.len(), self.moduli.len());

        match self.moduli.len() {
            0 => unreachable!(),
            1 => residues.to_vec(),
            _ => {
                let mut result = vec![T::ZERO; self.moduli_product.len()];
                let mut inter = vec![T::ZERO; self.moduli_product.len()];

                izip!(
                    residues,
                    &self.inv_punctured_moduli,
                    &self.punctured_moduli,
                    &self.moduli
                )
                .for_each(
                    |(&ri, &inv_mi, mi, &modulus): (&T, &ShoupFactor<T>, &Vec<T>, &M)| {
                        let product = inv_mi.factor_mul_modulo(ri, modulus.value_unchecked());
                        mi.slice_mul_value_inplace(product, &mut inter);
                        result.slice_add_modulo_assign(&inter, &self.moduli_product);
                    },
                );

                result
            }
        }
    }
}
