use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::{Rlev, RlevInfo};

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different powers
/// of a base, used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct NttRlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> NttRlev<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlev<S, T>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> NttRlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlev<S>`] that is initialized to zero.
    #[inline]
    pub fn zero(info: RlevInfo) -> Self {
        let len = info.decompose_length * info.poly_length * 2;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Creates a new [`NttRlev<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: ArrayBase::from_slice(converted_data),
        }
    }

    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &NttRlev<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &NttRlev<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
        self
    }

    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(mut self, ntt_table: &Table) -> Rlev<S>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();
        self.data.chunks_exact_mut(poly_length).for_each(|values| {
            ntt_table.inverse_transform_slice(values);
        });

        Rlev::new(self.data)
    }
}

impl<S, T> NttRlev<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`NttRlev<S>`] from bytes `data`.
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
    pub fn add_assign_element_wise<M, A>(&mut self, rhs: &NttRlev<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M, A>(&mut self, rhs: &NttRlev<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
    }
}

impl<S, T> NttRlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`NttRlev<S>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        converted_data.to_vec()
    }

    /// Converts [`NttRlev<S>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`NttRlev<S>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.byte_count()
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut NttRlev<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.add_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut NttRlev<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table, A>(&self, result: &mut Rlev<A>, ntt_table: &Table)
    where
        A: RawData<Elem = T> + DataMut,
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result
            .data
            .chunks_exact_mut(poly_length)
            .for_each(|values| {
                ntt_table.inverse_transform_slice(values);
            });
    }
}
