use primus_integer::{DataMut, RawData, UnsignedInteger, izip};
use primus_reduce::ops::{ReduceInv, ReduceMul};

use super::DcrtPolynomial;

#[inline]
unsafe fn invert_slice_unchecked<T, M>(values: &mut [T], modulus: M, scratch: &mut [T])
where
    T: UnsignedInteger,
    M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
{
    if values.is_empty() {
        return;
    }

    let prefix_products = &mut scratch[..values.len()];

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

impl<S, T> DcrtPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs the point-wise inverse for each modulus component.
    ///
    /// # Panics
    ///
    /// Panics if `scratch.len() < poly_length` or any value is zero.
    #[inline]
    pub fn inv<M>(mut self, poly_length: usize, moduli: &[M], scratch: &mut [T]) -> Self
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        self.inv_assign(poly_length, moduli, scratch);
        self
    }

    /// Performs the point-wise inverse for each modulus component in place.
    ///
    /// # Panics
    ///
    /// Panics if `scratch.len() < poly_length` or any value is zero.
    #[inline]
    pub fn inv_assign<M>(&mut self, poly_length: usize, moduli: &[M], scratch: &mut [T])
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        assert!(
            scratch.len() >= poly_length,
            "scratch length {} is smaller than polynomial length {}",
            scratch.len(),
            poly_length
        );
        assert!(
            self.0.iter().all(|v| !v.is_zero()),
            "point-wise inverse is undefined for zero values"
        );

        unsafe { self.inv_assign_unchecked(poly_length, moduli, scratch) }
    }

    /// Performs the point-wise inverse for each modulus component without
    /// checking whether any value is zero.
    ///
    /// # Safety
    ///
    /// - `scratch.len() >= poly_length`
    /// - every value in `self` is non-zero
    #[inline]
    pub unsafe fn inv_unchecked<M>(
        mut self,
        poly_length: usize,
        moduli: &[M],
        scratch: &mut [T],
    ) -> Self
    where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        unsafe { self.inv_assign_unchecked(poly_length, moduli, scratch) };
        self
    }

    /// Performs the point-wise inverse for each modulus component in place
    /// without checking whether any value is zero.
    ///
    /// Internally this uses batch inversion for each modulus block:
    /// one modular inverse plus `O(poly_length)` modular multiplications.
    ///
    /// # Safety
    ///
    /// - `scratch.len() >= poly_length`
    /// - every value in `self` is non-zero
    #[inline]
    pub unsafe fn inv_assign_unchecked<M>(
        &mut self,
        poly_length: usize,
        moduli: &[M],
        scratch: &mut [T],
    ) where
        M: Copy + ReduceMul<T, Output = T> + ReduceInv<T, Output = T>,
    {
        debug_assert!(scratch.len() >= poly_length);

        izip!(self.iter_each_modulus_mut(poly_length), moduli).for_each(|(poly, &modulus)| {
            unsafe { invert_slice_unchecked(poly, modulus, scratch) };
        });
    }
}
