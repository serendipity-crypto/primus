use integer::{BigIntegerOps, UnsignedInteger};
use num_traits::ConstOne;
use serde::{Deserialize, Serialize};

/// The basis for approximate signed decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct BigUintApproxSignedBasis<T: UnsignedInteger> {
    modulus: Vec<T>,
}

impl<T: UnsignedInteger> BigUintApproxSignedBasis<T> {
    #[inline]
    pub fn new(modulus: &[T], log_basis: u32, reverse_length: Option<usize>) -> Self {
        assert!(log_basis > 0 && T::BITS > log_basis);

        let basis = <T as ConstOne>::ONE << log_basis;
        let basis_minus_one = basis - <T as ConstOne>::ONE;
        let bits_count = modulus.slice_value_bits_count();

        todo!()
    }
}
