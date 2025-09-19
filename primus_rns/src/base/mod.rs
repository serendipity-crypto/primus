use integer::UnsignedInteger;
use itertools::Itertools;
use primus_factor::FactorMul;
use reduce::FieldAdapter;

use crate::RNSError;

pub struct RNSBase<T: UnsignedInteger, M: FieldAdapter<T>, F: FactorMul<T>> {
    pub moduli: Vec<M>,
    pub moduli_product: T,
    pub punctured_moduli: Vec<Vec<T>>,
    pub inv_punctured_moduli: Vec<F>,
}

impl<T: UnsignedInteger, M: FieldAdapter<T>, F: FactorMul<T>> RNSBase<T, M, F> {
    pub fn new(moduli: &[M]) -> Result<Self, RNSError> {
        let moduli_values = moduli
            .iter()
            .map(|m| m.value().unwrap())
            .collect::<Vec<_>>();

        if moduli_values
            .iter()
            .tuple_combinations()
            .any(|(&a, &b)| a.not_coprime(b))
        {
            return Err(RNSError::CoPrimeError);
        }
        todo!()
    }
}
