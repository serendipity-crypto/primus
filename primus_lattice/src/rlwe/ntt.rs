use primus_integer::{ByteCount, UnsignedInteger};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayMut, ArrayRef, NttPolynomialRef, PolynomialRef};
use primus_reduce::FieldContext;
use serde::{Deserialize, Serialize};

use crate::Rlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct NttRlwe<T: UnsignedInteger> {
    pub(crate) data: Vec<T>,
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Creates a new [`NttRlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: converted_data.to_vec(),
        }
    }

    /// Creates a new [`NttRlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.data.copy_from_slice(converted_data);
    }

    /// Converts [`NttRlwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_slice());

        converted_data.to_vec()
    }

    /// Converts [`NttRlwe<T>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_slice());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`NttRlwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Creates a new [`NttRlwe<T>`].
    #[inline]
    pub fn new(data: Vec<T>) -> Self {
        Self { data }
    }

    /// Creates a new [`NttRlwe<T>`] with reference of [`PolynomialRef<'_, T>`].
    #[inline]
    pub fn from_ref(a: PolynomialRef<'_, T>, b: PolynomialRef<'_, T>) -> Self {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: [a.0, b.0].concat(),
        }
    }

    /// Creates a new [`NttRlwe<T>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self {
            data: vec![T::ZERO; poly_length << 1],
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        ArrayMut(&mut self.data).set_zero();
    }

    /// Extracts mutable slice of `a` and `b` of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_mut_unchecked(mid) }
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(mut self, ntt_table: &Table) -> Rlwe<T>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a, b) = self.a_b_mut_slices();

        ntt_table.inverse_transform_slice(a);
        ntt_table.inverse_transform_slice(b);

        Rlwe::new(self.data)
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table>(&self, ntt_table: &Table, result: &mut Rlwe<T>)
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        result.data.copy_from_slice(&self.data);

        let (a, b) = result.a_b_mut_slices();

        ntt_table.inverse_transform_slice(a);
        ntt_table.inverse_transform_slice(b);
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: FieldContext<T>,
    {
        ArrayMut(&mut self.data).add_assign(ArrayRef(&rhs.data), modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: FieldContext<T>,
    {
        ArrayMut(&mut self.data).sub_assign(ArrayRef(&rhs.data), modulus);
        self
    }

    /// Performs an in-place element-wise modular addition `self += rhs`.
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        ArrayMut(&mut self.data).add_assign(ArrayRef(&rhs.data), modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        ArrayMut(&mut self.data).sub_assign(ArrayRef(&rhs.data), modulus);
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        ArrayRef(&self.data).add_inplace(ArrayRef(&rhs.data), ArrayMut(&mut result.data), modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        ArrayRef(&self.data).sub_inplace(ArrayRef(&rhs.data), ArrayMut(&mut result.data), modulus)
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<T>`] with another `polynomial` [`NttPolynomial<T>`].
    #[inline]
    pub fn mul_ntt_polynomial_assign<M>(&mut self, polynomial: NttPolynomialRef<'_, T>, modulus: M)
    where
        M: FieldContext<T>,
    {
        let (a, b) = self.a_b_mut_slices();
        // NttPolynomialMut(a).
        todo!()
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<T>`] with another `polynomial` [`NttPolynomial<T>`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M>(
        &self,
        polynomial: NttPolynomialRef<'_, T>,
        result: &mut Self,
        modulus: M,
    ) where
        M: FieldContext<T>,
    {
        // self.a.mul_inplace(polynomial, result.a_mut(), modulus);
        // self.b.mul_inplace(polynomial, result.b_mut(), modulus);

        todo!()
    }
}
