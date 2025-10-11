use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use crate::Rlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct NttRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T: UnsignedInteger> NttRlwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T: UnsignedInteger> NttRlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwe<S>`] with reference of [`NttPolynomial<A>`].
    #[inline]
    pub fn from_ref<A>(a: &NttPolynomial<A>, b: &NttPolynomial<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: ArrayBase::from_vec([a.0.as_ref(), b.0.as_ref()].concat()),
        }
    }

    /// Creates a new [`NttRlwe<S>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; poly_length << 1]),
        }
    }

    /// Creates a new [`NttRlwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: ArrayBase::from_slice(converted_data),
        }
    }

    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &NttRlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &NttRlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
        self
    }

    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(mut self, ntt_table: &Table) -> Rlwe<S>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a, b) = self.a_b_mut_slices();

        ntt_table.inverse_transform_slice(a);
        ntt_table.inverse_transform_slice(b);

        Rlwe::new(self.data)
    }
}

impl<S, T: UnsignedInteger> NttRlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.data.copy_from_slice(converted_data);
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.set_zero();
    }

    /// Extracts mutable slice of `a` and `b` of this [`NttRlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_mut_unchecked(mid) }
    }

    /// Performs an in-place element-wise modular addition `self += rhs`.
    #[inline]
    pub fn add_assign_element_wise<M, A>(&mut self, rhs: &NttRlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M, A>(&mut self, rhs: &NttRlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<S>`] with another `polynomial` [`NttPolynomial<A>`].
    #[inline]
    pub fn mul_ntt_polynomial_assign<M, A>(&mut self, polynomial: &NttPolynomial<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let (a, b) = self.a_b_mut_slices();

        NttPolynomial(ArrayBase(a)).mul_assign(polynomial, modulus);
        NttPolynomial(ArrayBase(b)).mul_assign(polynomial, modulus);
    }
}

impl<S, T: UnsignedInteger> NttRlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`NttRlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self) -> (&[T], &[T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Converts [`NttRlwe<S>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        converted_data.to_vec()
    }

    /// Converts [`NttRlwe<S>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`NttRlwe<S>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.byte_count()
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut NttRlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.add_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut NttRlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<S>`] with another `polynomial` [`NttPolynomial`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, A, B>(
        &self,
        polynomial: &NttPolynomial<A>,
        result: &mut NttRlwe<B>,
        modulus: M,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let (a0, b0) = self.a_b_slices();
        let (a1, b1) = result.a_b_mut_slices();

        NttPolynomial(ArrayBase(a0)).mul_inplace(
            polynomial,
            &mut NttPolynomial(ArrayBase(a1)),
            modulus,
        );
        NttPolynomial(ArrayBase(b0)).mul_inplace(
            polynomial,
            &mut NttPolynomial(ArrayBase(b1)),
            modulus,
        );
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table, A>(&self, result: &mut Rlwe<A>, ntt_table: &Table)
    where
        A: RawData<Elem = T> + DataMut,
        Table: NttTable<ValueT = T> + Ntt,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let (a, b) = result.a_b_mut_slices();

        ntt_table.inverse_transform_slice(a);
        ntt_table.inverse_transform_slice(b);
    }
}
