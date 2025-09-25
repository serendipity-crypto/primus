use integer::UnsignedInteger;
use num_traits::ConstOne;
use serde::{Deserialize, Serialize};

use super::{ScalarIter, SignedDecomposeIter};

/// The basis for approximate signed decomposition.
#[derive(Debug, Clone, Copy, Eq, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct ApproxSignedBasis<T: UnsignedInteger> {
    modulus: Option<T>,
    modulus_is_power_of_2: bool,
    basis: T,
    basis_minus_one: T,
    modulus_minus_basis: T,
    decompose_length: usize,
    value_bits: u32,
    init_carry_mask: Option<T>,
    carry_mask: T,
    log_basis: u32,
    drop_bits: u32,
    split_value: Option<T>,
    next_pow_of_2_sub_modulus: T,
}

impl<T: UnsignedInteger> PartialEq for ApproxSignedBasis<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.modulus == other.modulus
            && self.basis == other.basis
            && self.decompose_length == other.decompose_length
    }
}

impl<T: UnsignedInteger> ApproxSignedBasis<T> {
    #[inline]
    pub fn new(modulus: Option<T>, log_basis: u32, reverse_length: Option<usize>) -> Self {
        assert!(log_basis > 0);

        let basis = <T as ConstOne>::ONE << log_basis;
        let basis_minus_one = basis - <T as ConstOne>::ONE;

        let modulus_is_power_of_2;
        let value_bits;
        let modulus_minus_basis;

        if let Some(modulus) = modulus {
            if modulus.is_power_of_two() {
                modulus_is_power_of_2 = true;
                value_bits = modulus.trailing_zeros();
            } else {
                modulus_is_power_of_2 = false;
                value_bits = T::BITS - modulus.leading_zeros();
            }
            assert!(value_bits >= log_basis);
            modulus_minus_basis = modulus - basis;
        } else {
            assert!(T::BITS >= log_basis);
            modulus_is_power_of_2 = true;
            value_bits = T::BITS;
            modulus_minus_basis = T::MAX - basis_minus_one;
        }

        let decompose_length = value_bits / log_basis;
        let mut drop_bits = value_bits - decompose_length * log_basis;
        let mut decompose_length = decompose_length as usize;

        if let Some(reverse_len) = reverse_length {
            assert!(decompose_length >= reverse_len);
            decompose_length = reverse_len;
            drop_bits = value_bits - (reverse_len as u32) * log_basis;
        }

        assert!(decompose_length > 0);

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

        let mut split_value = None;
        let mut next_pow_of_2_sub_modulus = T::ZERO;
        if !modulus_is_power_of_2 {
            let modulus = modulus.unwrap();
            split_value = if log_basis == 1 {
                if drop_bits == 0 {
                    None
                } else {
                    let mut value = T::ZERO;
                    for _ in 0..decompose_length {
                        value <<= 1;
                        value |= T::ONE;
                    }
                    value <<= 1;
                    value |= T::ONE;
                    value <<= drop_bits - 1;
                    if value >= modulus { None } else { Some(value) }
                }
            } else {
                let mut value = T::ZERO;
                for _ in 0..decompose_length {
                    value <<= log_basis;
                    value |= basis_minus_one >> 1u32;
                }
                if drop_bits > 0 {
                    value <<= 1;
                    value |= T::ONE;
                    value <<= drop_bits - 1;
                } else {
                    value += T::ONE;
                }
                if value >= modulus { None } else { Some(value) }
            };

            next_pow_of_2_sub_modulus = (T::MAX >> (T::BITS - value_bits)) - (modulus - T::ONE);
        }

        Self {
            modulus,
            modulus_is_power_of_2,
            basis,
            basis_minus_one,
            modulus_minus_basis,
            value_bits,
            init_carry_mask,
            carry_mask,
            decompose_length,
            log_basis,
            drop_bits,
            split_value,
            next_pow_of_2_sub_modulus,
        }
    }

    /// Checks whether the modulus of this [`ApproxSignedBasis<T>`] is power of 2.
    #[inline]
    pub fn modulus_is_power_of_2(&self) -> bool {
        self.modulus_is_power_of_2
    }

    /// Returns the value bits of values in `[0, modulus - 1]`.
    #[inline]
    pub fn value_bits(&self) -> u32 {
        self.value_bits
    }

    /// Returns the decompose length of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_length(&self) -> usize {
        self.decompose_length
    }

    /// Returns the basis value of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn basis_value(&self) -> T {
        self.basis
    }

    /// Returns the basis minus one of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn basis_minus_one(&self) -> T {
        self.basis_minus_one
    }

    /// Returns the log basis of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn log_basis(&self) -> u32 {
        self.log_basis
    }

    /// Returns the drop bits of this [`ApproxSignedBasis<T>`].
    ///
    /// This means some bits of the value will be dropped
    /// according to approximate signed decomposition.
    #[inline]
    pub fn drop_bits(&self) -> u32 {
        self.drop_bits
    }

    /// Returns the init carry mask of this [`ApproxSignedBasis<T>`].
    ///
    /// This value is used for generating the initial carry for decomposition.
    #[inline]
    pub fn init_carry_mask(&self) -> Option<T> {
        self.init_carry_mask
    }

    /// Returns an iterator over the signed decomposition operators of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_iter(&self) -> SignedDecomposeIter<T> {
        SignedDecomposeIter::<T> {
            length: self.decompose_length,
            value_chunk_mask: self.basis_minus_one << self.drop_bits,
            mask_shl_bits: self.log_basis,
            value_shr_bits: self.drop_bits,
            carry_mask: self.carry_mask,
            basis_minus_one: self.basis_minus_one,
            modulus_minus_basis: self.modulus_minus_basis,
        }
    }

    /// Returns an iterator over scalars of this [`ApproxSignedBasis<T>`].
    #[inline]
    pub fn scalar_iter(&self) -> ScalarIter<T> {
        ScalarIter::new(
            T::ONE << self.drop_bits,
            self.decompose_length,
            self.log_basis,
        )
    }

    /// Init carry and adjusted value for a value.
    #[inline]
    pub fn init_value_carry(&self, value: T) -> (T, bool) {
        let mut adjust = value;
        if let Some(split) = self.split_value {
            if value >= split {
                adjust += self.next_pow_of_2_sub_modulus;
            }
        }

        (
            adjust,
            match self.init_carry_mask {
                Some(mask) => !(adjust & mask).is_zero(),
                None => false,
            },
        )
    }

    /// Init carries and adjusted values for a slice and store the adjusted values back to `values`.
    #[inline]
    pub fn init_value_carry_slice_inplace(&self, values: &mut [T], carries: &mut [bool]) {
        if let Some(split) = self.split_value {
            values.iter_mut().for_each(|value| {
                if *value >= split {
                    *value += self.next_pow_of_2_sub_modulus;
                }
            })
        }

        match self.init_carry_mask {
            Some(mask) => values.iter().zip(carries).for_each(|(&value, carry)| {
                *carry = !(value & mask).is_zero();
            }),
            None => carries.fill(false),
        };
    }

    /// Init carries and adjusted values for a slice.
    #[inline]
    pub fn init_value_carry_slice(
        &self,
        values: &[T],
        adjust_values: &mut [T],
        carries: &mut [bool],
    ) {
        if let Some(split) = self.split_value {
            adjust_values
                .iter_mut()
                .zip(values)
                .for_each(|(adjust_value, &value)| {
                    if value >= split {
                        *adjust_value = value + self.next_pow_of_2_sub_modulus;
                    } else {
                        *adjust_value = value;
                    }
                })
        } else {
            adjust_values.copy_from_slice(values);
        }

        match self.init_carry_mask {
            Some(mask) => adjust_values
                .iter()
                .zip(carries)
                .for_each(|(&value, carry)| {
                    *carry = !(value & mask).is_zero();
                }),
            None => carries.fill(false),
        };
    }
}
