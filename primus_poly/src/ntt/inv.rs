use primus_integer::{DataMut, RawData, UnsignedInteger};
use primus_reduce::ops::{ReduceInv, ReduceMul};

use super::NttPolynomial;

impl<S, T> NttPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the point-wise inverse in the NTT domain.
    ///
    /// # Panics
    ///
    /// Panics if `scratch.len() < self.poly_length()` or any value is zero.
    #[inline]
    pub fn inv<M>(mut self, modulus: M, scratch: &mut [T]) -> Self
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        self.inv_assign(modulus, scratch);
        self
    }

    /// Performs the point-wise inverse in the NTT domain in place.
    ///
    /// # Panics
    ///
    /// Panics if `scratch.len() < self.poly_length()` or any value is zero.
    #[inline]
    pub fn inv_assign<M>(&mut self, modulus: M, scratch: &mut [T])
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        assert!(
            scratch.len() >= self.poly_length(),
            "scratch length {} is smaller than polynomial length {}",
            scratch.len(),
            self.poly_length()
        );
        assert!(
            self.iter().all(|v| !v.is_zero()),
            "point-wise inverse is undefined for zero values"
        );

        unsafe { self.inv_assign_unchecked(modulus, scratch) }
    }

    /// Performs the point-wise inverse in the NTT domain without checking
    /// whether any value is zero.
    ///
    /// # Safety
    ///
    /// - `scratch.len() >= self.poly_length()`
    /// - every value in `self` is non-zero
    #[inline]
    pub unsafe fn inv_unchecked<M>(mut self, modulus: M, scratch: &mut [T]) -> Self
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        unsafe { self.inv_assign_unchecked(modulus, scratch) };
        self
    }

    /// Performs the point-wise inverse in the NTT domain in place without
    /// checking whether any value is zero.
    ///
    /// Internally this uses batch inversion:
    /// one modular inverse plus `O(n)` modular multiplications.
    ///
    /// # Safety
    ///
    /// - `scratch.len() >= self.poly_length()`
    /// - every value in `self` is non-zero
    #[inline]
    pub unsafe fn inv_assign_unchecked<M>(&mut self, modulus: M, scratch: &mut [T])
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        let poly_length = self.poly_length();
        debug_assert!(scratch.len() >= poly_length);

        if poly_length == 0 {
            return;
        }

        let values = self.as_mut_slice();
        let prefix_products = &mut scratch[..poly_length];

        let mut total_product = T::ONE;
        for (prefix_product, &value) in prefix_products.iter_mut().zip(values.iter()) {
            *prefix_product = total_product;
            total_product = modulus.reduce_mul(total_product, value);
        }

        let mut suffix_inverse = modulus.reduce_inv(total_product);

        for (value, prefix_product) in values
            .iter_mut()
            .rev()
            .zip(prefix_products.iter().rev().copied())
        {
            let current_value = *value;
            *value = modulus.reduce_mul(prefix_product, suffix_inverse);
            suffix_inverse = modulus.reduce_mul(suffix_inverse, current_value);
        }
    }
}
