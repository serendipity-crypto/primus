use integer::{BigIntegerOps, UnsignedInteger};
use num_traits::ConstOne;
use serde::{Deserialize, Serialize};

use super::{BigUintScalarIter, BigUintSignedDecomposeIter, ValueMask};

/// The basis for approximate signed decomposition.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct BigUintApproxSignedBasis<T: UnsignedInteger> {
    modulus: Vec<T>,
    basis: T,
    basis_minus_one: T,
    decompose_length: usize,
    log_basis: u32,
    drop_bits: u32,
    init_carry_mask: Option<(usize, T)>,
    init_value_mask: ValueMask<T>,
    carry_mask: T,
    split_value: Option<Vec<T>>,
    modulus_minus_basis: Vec<T>,
    next_pow_of_2_sub_modulus: Vec<T>,
}

impl<T: UnsignedInteger> BigUintApproxSignedBasis<T> {
    #[inline]
    pub fn new(modulus: &[T], log_basis: u32, reverse_length: Option<usize>) -> Self {
        // FIXME: log_basis may be bigger than T::BITS
        assert!(!modulus.last().unwrap().is_zero());
        assert!(log_basis > 0 && T::BITS > log_basis);

        let len = modulus.len();
        let unused_bits = modulus[len - 1].leading_zeros();

        let basis = <T as ConstOne>::ONE << log_basis;
        let basis_minus_one = basis - <T as ConstOne>::ONE;
        let modulus_bits_count = T::BITS * (len as u32) - unused_bits;
        let decompose_length = modulus_bits_count / log_basis;
        let mut drop_bits = modulus_bits_count - decompose_length * log_basis;
        let mut decompose_length = decompose_length as usize;

        if let Some(reverse_len) = reverse_length {
            assert!(decompose_length >= reverse_len);
            decompose_length = reverse_len;
            drop_bits = modulus_bits_count - (reverse_len as u32) * log_basis;
        }

        assert!(decompose_length > 0);

        let init_carry_mask = if drop_bits > 0 {
            let bits = drop_bits - 1;
            let (index, bits) = (bits / T::BITS, bits % T::BITS);
            Some((index as usize, T::ONE << bits))
        } else {
            None
        };

        let init_value_mask = ValueMask::new(basis_minus_one, drop_bits);

        let carry_mask = if log_basis == 1 {
            T::ONE << 1u32
        } else {
            (T::ONE << log_basis) | (T::ONE << (log_basis - 1))
        };

        let split_value: Option<Vec<T>> = if log_basis == 1 {
            if drop_bits == 0 {
                None
            } else {
                let mut value = vec![T::ZERO; modulus.len()];
                for _ in 0..decompose_length {
                    value.slice_left_shift_assign(1);
                    value[0] |= T::ONE;
                }
                value.slice_left_shift_assign(1);
                value[0] |= T::ONE;
                value.slice_left_shift_assign(drop_bits - 1);
                if value.slice_cmp(modulus).is_ge() {
                    None
                } else {
                    Some(value)
                }
            }
        } else {
            let mut value = vec![T::ZERO; modulus.len()];
            for _ in 0..decompose_length {
                value.slice_left_shift_assign(log_basis);
                value[0] |= basis_minus_one >> 1u32;
            }
            if drop_bits > 0 {
                value.slice_left_shift_assign(1);
                value[0] |= T::ONE;
                value.slice_left_shift_assign(drop_bits - 1);
            } else {
                let carry = value.slice_add_value_assign(T::ONE);
                assert!(!carry);
            }

            if value.slice_cmp(modulus).is_ge() {
                None
            } else {
                Some(value)
            }
        };

        let mut modulus_minus_basis = modulus.to_vec();
        let borrow = modulus_minus_basis.slice_sub_value_assign(basis);
        assert!(!borrow);

        let next_pow_of_2_sub_modulus: Vec<T> = {
            let mut next_pow_of_2_minus_one = vec![T::MAX; modulus.len()];
            next_pow_of_2_minus_one[len - 1] >>= unused_bits;

            let mut modulus_minus_one = modulus.to_vec();
            let _ = modulus_minus_one.slice_sub_value_assign(T::ONE);

            let borrow = next_pow_of_2_minus_one.slice_sub_assign(&modulus_minus_one);
            assert!(!borrow);
            next_pow_of_2_minus_one
        };

        Self {
            modulus: modulus.to_vec(),
            basis,
            basis_minus_one,
            decompose_length,
            log_basis,
            drop_bits,
            init_carry_mask,
            init_value_mask,
            carry_mask,
            split_value,
            modulus_minus_basis,
            next_pow_of_2_sub_modulus,
        }
    }

    /// Returns a reference to the modulus of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn modulus(&self) -> &[T] {
        &self.modulus
    }

    /// Returns the basis of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn basis_value(&self) -> T {
        self.basis
    }

    /// Returns the basis minus one of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn basis_minus_one(&self) -> T {
        self.basis_minus_one
    }

    /// Returns the decompose length of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_length(&self) -> usize {
        self.decompose_length
    }

    /// Returns the log basis of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn log_basis(&self) -> u32 {
        self.log_basis
    }

    /// Returns the drop bits of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn drop_bits(&self) -> u32 {
        self.drop_bits
    }

    /// Returns the init carry mask of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn init_carry_mask(&self) -> Option<(usize, T)> {
        self.init_carry_mask
    }

    /// Returns the carry mask of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn carry_mask(&self) -> T {
        self.carry_mask
    }

    /// Returns the split value of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn split_value(&self) -> Option<&Vec<T>> {
        self.split_value.as_ref()
    }

    /// Returns a reference to the modulus minus basis of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn modulus_minus_basis(&self) -> &[T] {
        &self.modulus_minus_basis
    }

    /// Returns a reference to the next pow of 2 sub modulus of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn next_pow_of_2_sub_modulus(&self) -> &[T] {
        &self.next_pow_of_2_sub_modulus
    }

    /// Returns an iterator over the signed decomposition operators of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn decompose_iter(&self) -> BigUintSignedDecomposeIter<T> {
        BigUintSignedDecomposeIter::<T> {
            length: self.decompose_length,
            value_mask: self.init_value_mask,
            mask_shl_bits: self.log_basis,
            carry_mask: self.carry_mask,
            basis_minus_one: self.basis_minus_one,
            modulus_minus_basis: self.modulus_minus_basis.clone(),
        }
    }

    /// Returns an iterator over scalars of this [`BigUintApproxSignedBasis<T>`].
    #[inline]
    pub fn scalar_iter(&self) -> BigUintScalarIter<T> {
        let mut scalar = vec![T::ZERO; self.modulus.len()];
        scalar[0] = T::ONE;
        scalar.slice_left_shift_assign(self.drop_bits);
        BigUintScalarIter::new(scalar, self.decompose_length, self.log_basis)
    }

    /// Init carry and adjusted value for a value.
    #[inline]
    pub fn init_value_carry(&self, value: &[T]) -> (Vec<T>, bool) {
        let mut adjust = value.to_vec();
        if let Some(split) = &self.split_value {
            if value.slice_cmp(&split).is_ge() {
                let _ = adjust.slice_add_assign(&self.next_pow_of_2_sub_modulus);
            }
        }

        let carry = match self.init_carry_mask {
            Some((i, mask)) => !(adjust[i] & mask).is_zero(),
            None => false,
        };

        (adjust, carry)
    }

    /// Init carries and adjusted values for a slice and store the adjusted values back to `values`.
    #[inline]
    pub fn init_value_carry_slice_inplace(
        &self,
        values: &mut [T],
        carries: &mut [bool],
        value_chunk_size: usize,
    ) {
        if let Some(split) = &self.split_value {
            values.chunks_exact_mut(value_chunk_size).for_each(|value| {
                if value.slice_cmp(&split).is_ge() {
                    let _ = value.slice_add_assign(&self.next_pow_of_2_sub_modulus);
                }
            })
        }

        match self.init_carry_mask {
            Some((i, mask)) => values
                .chunks_exact_mut(value_chunk_size)
                .zip(carries)
                .for_each(|(value, carry)| {
                    *carry = !(value[i] & mask).is_zero();
                }),
            None => carries.fill(false),
        };
    }

    /// Init carries and adjusted values for a slice.
    #[inline]
    pub fn init_value_carry_slice(
        &self,
        values: &[T],
        carries: &mut [bool],
        adjust_values: &mut [T],
        value_chunk_size: usize,
    ) {
        if let Some(split) = &self.split_value {
            adjust_values
                .chunks_exact_mut(value_chunk_size)
                .zip(values.chunks_exact(value_chunk_size))
                .for_each(|(adjust_value, value)| {
                    adjust_value.copy_from_slice(value);
                    if value.slice_cmp(&split).is_ge() {
                        let _ = adjust_value.slice_add_assign(&self.next_pow_of_2_sub_modulus);
                    }
                })
        }

        match self.init_carry_mask {
            Some((i, mask)) => adjust_values
                .chunks_exact_mut(value_chunk_size)
                .zip(carries)
                .for_each(|(value, carry)| {
                    *carry = !(value[i] & mask).is_zero();
                }),
            None => carries.fill(false),
        };
    }
}
