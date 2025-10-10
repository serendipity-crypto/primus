use primus_integer::{ByteCount, UnsignedInteger};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayMut, ArrayRef, PolynomialRef};
use primus_reduce::FieldContext;
use serde::{Deserialize, Serialize};

use crate::NttRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct Rlwe<T: UnsignedInteger> {
    pub(crate) data: Vec<T>,
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Creates a new [`Rlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: converted_data.to_vec(),
        }
    }

    /// Creates a new [`Rlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.data.copy_from_slice(converted_data);
    }

    /// Converts [`Rlwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_slice());

        converted_data.to_vec()
    }

    /// Converts [`Rlwe<T>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_slice());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`Rlwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.len() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Creates a new [`Rlwe<T>`].
    #[inline]
    pub fn new(data: Vec<T>) -> Self {
        Self { data }
    }

    /// Creates a new [`Rlwe<T>`] with reference of [`Polynomial<T>`].
    #[inline]
    pub fn from_ref(a: PolynomialRef<'_, T>, b: PolynomialRef<'_, T>) -> Self {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: [a.0, b.0].concat(),
        }
    }

    /// Creates a new [`Rlwe<T>`] that is initialized to zero,
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

    /// Extracts mutable slice of `a` and `b` of this [`Rlwe<T>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_mut_unchecked(mid) }
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(mut self, ntt_table: &Table) -> NttRlwe<T>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a, b) = self.a_b_mut_slices();

        ntt_table.transform_slice(a);
        ntt_table.transform_slice(b);

        NttRlwe::new(self.data)
    }

    /// ntt transform
    #[inline]
    pub fn to_ntt_form_inplace<Table>(&self, ntt_table: &Table, result: &mut NttRlwe<T>)
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        result.data.copy_from_slice(&self.data);

        let (a, b) = result.a_b_mut_slices();

        ntt_table.transform_slice(a);
        ntt_table.transform_slice(b);
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
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
}

impl<T: UnsignedInteger> Rlwe<T> {}
