use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use crate::{GlweInfo, NttGlwe};

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct Glwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> Glwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`Glwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> Glwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Glwe<S>`] that is initialized to zero.
    #[inline]
    pub fn zero(info: GlweInfo) -> Self {
        let len = (info.dimension + 1) * info.poly_length;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Creates a new [`Glwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: ArrayBase::from_slice(converted_data),
        }
    }

    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &Glwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &Glwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
        self
    }

    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(mut self, ntt_table: &Table) -> NttGlwe<S>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();
        self.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
        });

        NttGlwe::new(self.data)
    }
}

impl<S, T> Glwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`Glwe<S>`] from bytes `data`.
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

    /// Performs an in-place element-wise modular addition `self += rhs`.
    #[inline]
    pub fn add_assign_element_wise<M, A>(&mut self, rhs: &Glwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M, A>(&mut self, rhs: &Glwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
    }
}

impl<S, T> Glwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`Glwe<S>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        converted_data.to_vec()
    }

    /// Converts [`Glwe<S>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`Glwe<S>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.byte_count()
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut Glwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.add_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut Glwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs a multiplication on the `self` [`Glwe<S>`] with another `ntt_polynomial` [`NttPolynomial<A>`],
    /// store the result into `result` [`NttGlwe<B>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table, A, B>(
        &self,
        ntt_polynomial: &NttPolynomial<A>,
        result: &mut NttGlwe<B>,
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
    pub fn to_ntt_form_inplace<Table, A>(&self, result: &mut NttGlwe<A>, ntt_table: &Table)
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
