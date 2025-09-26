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
        // FIXME: log_basis may be bigger than T::BITS
        assert!(log_basis > 0 && T::BITS > log_basis);

        let basis = <T as ConstOne>::ONE << log_basis;
        let basis_minus_one = basis - <T as ConstOne>::ONE;
        let value_bits_count = modulus.slice_value_bits_count();
        let decompose_length = value_bits_count / log_basis;
        let mut drop_bits = value_bits_count - decompose_length * log_basis;
        let mut decompose_length = decompose_length as usize;

        if let Some(reverse_len) = reverse_length {
            assert!(decompose_length >= reverse_len);
            decompose_length = reverse_len;
            drop_bits = value_bits_count - (reverse_len as u32) * log_basis;
        }

        assert!(decompose_length > 0);

        // FIXME: drop_bits may be bigger than T::BITS
        let init_carry_mask = if drop_bits > 0 {
            Some(<T as ConstOne>::ONE << (drop_bits - 1))
        } else {
            None
        };

        let carry_mask = if log_basis == 1 {
            T::ONE << 1u32
        } else {
            (T::ONE << log_basis) | (T::ONE << (log_basis - 1))
        };

        let split_value: Option<Vec<T>> = if log_basis == 1 {
            if drop_bits == 0 { None } else { todo!() }
        } else {
            todo!()
        };

        let next_pow_of_2_sub_modulus: Vec<T> = todo!();

        todo!()
    }
}
