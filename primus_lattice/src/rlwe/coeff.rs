use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, Polynomial, RawData};
use primus_reduce::FieldContext;

use crate::NttRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct Rlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`Rlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Rlwe<S>`] with reference of [`Polynomial<A>`].
    #[inline]
    pub fn from_ref<A>(a: &Polynomial<A>, b: &Polynomial<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: ArrayBase::from_vec([a.0.as_ref(), b.0.as_ref()].concat()),
        }
    }

    /// Creates a new [`Rlwe<S>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; poly_length << 1]),
        }
    }

    /// Creates a new [`Rlwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: ArrayBase::from_slice(converted_data),
        }
    }

    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &Rlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &Rlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
        self
    }

    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(mut self, ntt_table: &Table) -> NttRlwe<S>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();
        self.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
        });

        NttRlwe::new(self.data)
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`Rlwe<S>`] from bytes `data`.
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

    /// Extracts mutable slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.0.split_at_mut_unchecked(mid) }
    }

    /// Performs an in-place element-wise modular addition `self += rhs`.
    #[inline]
    pub fn add_assign_element_wise<M, A>(&mut self, rhs: &Rlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M, A>(&mut self, rhs: &Rlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self) -> (&[T], &[T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Converts [`Rlwe<S>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        converted_data.to_vec()
    }

    /// Converts [`Rlwe<S>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`Rlwe<S>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.byte_count()
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut Rlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.add_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut Rlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs a multiplication on the `self` [`Rlwe<S>`] with another `ntt_polynomial` [`NttPolynomial<A>`],
    /// store the result into `result` [`NttRlwe<B>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table, A, B>(
        &self,
        ntt_polynomial: &NttPolynomial<A>,
        result: &mut NttRlwe<B>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = ntt_table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
            NttPolynomial(ArrayBase(poly)).mul_assign(ntt_polynomial, modulus);
        });
    }

    /// ntt transform
    #[inline]
    pub fn to_ntt_form_inplace<Table, A>(&self, result: &mut NttRlwe<A>, ntt_table: &Table)
    where
        A: RawData<Elem = T> + DataMut,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
        });
    }
}
